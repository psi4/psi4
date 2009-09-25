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
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);

    /*** alpha-alpha-beta term 1 ***/ 
    /* quantities to compute X3 */
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_init(&WAbEi, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
    dpd_buf4_init(&WMbIj, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
    /* quantities to compute sigma */
    dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&WmAEf, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAmEf (mA,Ef)");
    dpd_buf4_init(&WMnIe, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
  
         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    if (params.t3_Ws_incore)
      cc3_sigma_RHF_ic(&CMnEf, &WAbEi, &WMbIj, 1,  &Dints, &SIA, 
        1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
        moinfo.virtpi, moinfo.vir_off, omega, outfile, params.nthreads,
        params.newtrips);
    else
      cc3_sigma_RHF(&CMnEf, &WAbEi, &WMbIj, 1,  &Dints, &SIA, 
        1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
        moinfo.virtpi, moinfo.vir_off, omega, outfile, params.newtrips);
  
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&WMbIj);
    dpd_buf4_close(&Dints);
    dpd_file2_close(&FME);
    dpd_buf4_close(&WmAEf);
    dpd_buf4_close(&WMnIe);

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
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&WAbEi, CC3_HC1ET1, C_irr, 10, 5, 10, 5, 0, "Ht_WAbEi (iE,bA)");
    dpd_buf4_init(&WMbIj, CC3_HC1ET1, C_irr, 0, 10, 0, 10, 0, "Ht_WMbIj (Ij,Mb)");
  
    dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&WmAEf, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAmEf (mA,Ef)");
    dpd_buf4_init(&WMnIe, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
  
           /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
              * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */
  
    if (params.t3_Ws_incore)
      cc3_sigma_RHF_ic(&tIjAb, &WAbEi, &WMbIj, 1,  &Dints, &SIA, 
        1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
        moinfo.virtpi, moinfo.vir_off, omega, outfile, params.nthreads,
        params.newtrips);
    else
      cc3_sigma_RHF(&tIjAb, &WAbEi, &WMbIj, 1,  &Dints, &SIA,
         1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
         moinfo.virtpi, moinfo.vir_off, omega, outfile, params.newtrips);
  
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&WMbIj);
    dpd_buf4_close(&Dints);
    dpd_file2_close(&FME);
    dpd_buf4_close(&WmAEf);
    dpd_buf4_close(&WMnIe);
  
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
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&WAbEi, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
    dpd_buf4_init(&WMbIj, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
  
    dpd_file2_init(&FME, CC3_HC1, C_irr, 0, 1, "HC1 FME");
    dpd_buf4_init(&WmAEf, CC3_HC1, C_irr, 10, 5, 10, 5, 0, "HC1 WAmEf (mA,Ef)");
    dpd_buf4_init(&WMnIe, CC3_HC1, C_irr, 0, 10, 0, 10, 0, "HC1 WMnIe (Mn,Ie)");
  
           /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */
  
    if (params.t3_Ws_incore)
      cc3_sigma_RHF_ic(&tIjAb, &WAbEi, &WMbIj, 0, NULL, NULL, 
        1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
        moinfo.virtpi, moinfo.vir_off, 0.0, outfile, params.nthreads,
        params.newtrips);
    else
    cc3_sigma_RHF(&tIjAb, &WAbEi, &WMbIj, 0, NULL, NULL,
       1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
       moinfo.virtpi, moinfo.vir_off, 0.0, outfile, params.newtrips);
  
    dpd_buf4_close(&tIjAb); 
    dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&WMbIj);
    dpd_file2_close(&FME);
    dpd_buf4_close(&WmAEf);
    dpd_buf4_close(&WMnIe);

#ifdef EOM_DEBUG
    dpd_file2_close(&SIA);
    dpd_buf4_close(&SIjAb);
    check_sum("<Psi|H'<T|(Uhat T2)c|0>|T>/(w-wt)", i, C_irr);
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
#endif

    dpd_file2_close(&SIA);
    dpd_buf4_close(&SIjAb);
  }
  else if (params.eom_ref == 1) { /* ROHF */
  }
  else if (params.eom_ref == 2) { /* UHF */
    /* open all sigma (output) files */
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

    /*** alpha-alpha-alpha term 1 */

    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);
    dpd_buf4_init(&WABEI, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    dpd_buf4_init(&WMBIJ, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    dpd_buf4_init(&DIJAB_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&WAMEF, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    dpd_buf4_init(&WMNIE, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");

         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    cc3_sigma_UHF_AAA(&CMNEF, &WABEI, &WMBIJ, 1, &DIJAB_anti, &SIA,
        1, &FME, &WAMEF, &WMNIE, &SIJAB, moinfo.aoccpi, moinfo.aocc_off,
        moinfo.avirtpi, moinfo.avir_off, omega, outfile);

    dpd_buf4_close(&CMNEF);
    dpd_buf4_close(&WABEI);
    dpd_buf4_close(&WMBIJ);
    dpd_buf4_close(&DIJAB_anti);
    dpd_file2_close(&FME);
    dpd_buf4_close(&WAMEF);
    dpd_buf4_close(&WMNIE);

    /*** beta-beta-beta term 1 */

    sprintf(lbl, "%s %d", "Cmnef", i);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 10, 15, 12, 17, 0, lbl);
    dpd_buf4_init(&Wabei, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    dpd_buf4_init(&Wmbij, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    dpd_buf4_init(&Dijab_anti, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&Wamef, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    dpd_buf4_init(&Wmnie, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");

         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    cc3_sigma_UHF_BBB(&Cmnef, &Wabei, &Wmbij, 1, &Dijab_anti, &Sia,
        1, &Fme, &Wamef, &Wmnie, &Sijab, moinfo.boccpi, moinfo.bocc_off,
        moinfo.bvirtpi, moinfo.bvir_off, omega, outfile);

    dpd_buf4_close(&Cmnef);
    dpd_buf4_close(&Wabei);
    dpd_buf4_close(&Wmbij);
    dpd_buf4_close(&Dijab_anti);
    dpd_file2_close(&Fme);
    dpd_buf4_close(&Wamef);
    dpd_buf4_close(&Wmnie);

    /*** alpha-alpha-beta term 1 */ 

    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    dpd_buf4_init(&CmNeF, EOM_TMP, C_irr, 23, 29, 23, 29, 0, "CmNeF");

    dpd_buf4_init(&WABEI, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    dpd_buf4_init(&WaBeI, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
    dpd_buf4_init(&WAbEi, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
    dpd_buf4_init(&WMBIJ, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    dpd_buf4_init(&WMbIj, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
    dpd_buf4_init(&WmBiJ, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

    dpd_buf4_init(&DIJAB_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    dpd_buf4_init(&DIjAb, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&WAMEF, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    dpd_buf4_init(&WaMeF, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");
    dpd_buf4_init(&WAmEf, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");
    dpd_buf4_init(&WMNIE, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");
    dpd_buf4_init(&WMnIe, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_init(&WmNiE, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");

         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    cc3_sigma_UHF_AAB(&CMNEF, &CMnEf, &CmNeF, &WABEI, &WaBeI, &WAbEi,
       &WMBIJ, &WMbIj, &WmBiJ, 1,  &DIJAB_anti, &DIjAb, &SIA, &Sia,
       1, &FME, &Fme, &WAMEF, &WaMeF, &WAmEf, &WMNIE, &WMnIe, &WmNiE,
       &SIJAB, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
       moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
       moinfo.bvir_off, omega, outfile);

    dpd_buf4_close(&CMNEF); dpd_buf4_close(&CMnEf); dpd_buf4_close(&CmNeF);
    dpd_buf4_close(&WABEI); dpd_buf4_close(&WaBeI); dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&WMBIJ); dpd_buf4_close(&WMbIj); dpd_buf4_close(&WmBiJ);
    dpd_buf4_close(&DIJAB_anti); dpd_buf4_close(&DIjAb);
    dpd_file2_close(&FME); dpd_file2_close(&Fme);
    dpd_buf4_close(&WAMEF); dpd_buf4_close(&WaMeF); dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&WMNIE); dpd_buf4_close(&WMnIe); dpd_buf4_close(&WmNiE);

    /*** beta-beta-alpha term 1 */ 

    sprintf(lbl, "%s %d", "Cmnef", i);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 10, 15, 12, 17, 0, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    dpd_buf4_init(&CmNeF, EOM_TMP, C_irr, 23, 29, 23, 29, 0, "CmNeF");

    dpd_buf4_init(&Wabei, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    dpd_buf4_init(&WaBeI, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
    dpd_buf4_init(&WAbEi, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
    dpd_buf4_init(&Wmbij, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    dpd_buf4_init(&WMbIj, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
    dpd_buf4_init(&WmBiJ, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

    dpd_buf4_init(&Dijab_anti, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    dpd_buf4_init(&DiJaB, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&Wamef, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    dpd_buf4_init(&WaMeF, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");
    dpd_buf4_init(&WAmEf, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");
    dpd_buf4_init(&Wmnie, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");
    dpd_buf4_init(&WMnIe, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_init(&WmNiE, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");

         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    cc3_sigma_UHF_BBA(&Cmnef, &CMnEf, &CmNeF, &Wabei, &WaBeI, &WAbEi,
      &Wmbij, &WMbIj, &WmBiJ, 1, &Dijab_anti, &DiJaB, &SIA, &Sia,
      1, &FME, &Fme, &Wamef, &WaMeF, &WAmEf, &Wmnie, &WMnIe, &WmNiE,
      &Sijab, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
      moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
      moinfo.bvir_off, omega, outfile);

    dpd_buf4_close(&Cmnef); dpd_buf4_close(&CMnEf); dpd_buf4_close(&CmNeF);
    dpd_buf4_close(&Wabei); dpd_buf4_close(&WaBeI); dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&Wmbij); dpd_buf4_close(&WMbIj); dpd_buf4_close(&WmBiJ);
    dpd_buf4_close(&Dijab_anti); dpd_buf4_close(&DiJaB);
    dpd_file2_close(&FME); dpd_file2_close(&Fme);
    dpd_buf4_close(&Wamef); dpd_buf4_close(&WaMeF); dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&Wmnie); dpd_buf4_close(&WMnIe); dpd_buf4_close(&WmNiE);

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

    dpd_buf4_init(&TIJAB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&WABEI, CC3_HC1ET1, C_irr, 20, 5, 20, 7, 0, "Ht_WABEI (IE,B>A)");
    dpd_buf4_init(&WMBIJ, CC3_HC1ET1, C_irr, 0, 20, 2, 20, 0, "Ht_WMBIJ (I>J,MB)");
    dpd_buf4_init(&DIJAB_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&WAMEF, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    dpd_buf4_init(&WMNIE, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");

         /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */

    cc3_sigma_UHF_AAA(&TIJAB, &WABEI, &WMBIJ, 1, &DIJAB_anti, &SIA,
        1, &FME, &WAMEF, &WMNIE, &SIJAB, moinfo.aoccpi, moinfo.aocc_off, 
        moinfo.avirtpi, moinfo.avir_off, omega, outfile);

    dpd_buf4_close(&TIJAB);
    dpd_buf4_close(&WABEI);
    dpd_buf4_close(&WMBIJ);
    dpd_buf4_close(&DIJAB_anti);
    dpd_file2_close(&FME);
    dpd_buf4_close(&WAMEF);
    dpd_buf4_close(&WMNIE);

    /*** beta-beta-beta term 2 */

    dpd_buf4_init(&Tijab, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&Wabei, CC3_HC1ET1, C_irr, 30, 15, 30, 17, 0, "Ht_Wabei (ie,b>a)");
    dpd_buf4_init(&Wmbij, CC3_HC1ET1, C_irr, 10, 30, 12, 30, 0, "Ht_Wmbij (i>j,mb)");
    dpd_buf4_init(&Dijab_anti, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&Wamef, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    dpd_buf4_init(&Wmnie, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");

         /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */

    cc3_sigma_UHF_BBB(&Tijab, &Wabei, &Wmbij, 1, &Dijab_anti, &Sia,
        1, &Fme, &Wamef, &Wmnie, &Sijab, moinfo.boccpi, moinfo.bocc_off,
        moinfo.bvirtpi, moinfo.bvir_off, omega, outfile);

    dpd_buf4_close(&Tijab);
    dpd_buf4_close(&Wabei);
    dpd_buf4_close(&Wmbij);
    dpd_buf4_close(&Dijab_anti);
    dpd_file2_close(&Fme);
    dpd_buf4_close(&Wamef);
    dpd_buf4_close(&Wmnie);

    /*** do alpha-alpha-beta term 2 */

    dpd_buf4_init(&TIJAB, CC_TAMPS, 0,  0,  5,  2,  7, 0, "tIJAB");
    dpd_buf4_init(&TIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&TiJaB, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

    dpd_buf4_init(&WABEI, CC3_HC1ET1, C_irr, 20, 5, 20, 7, 0, "Ht_WABEI (IE,B>A)");
    dpd_buf4_init(&WaBeI, CC3_HC1ET1, C_irr, 24, 28, 24, 28, 0, "Ht_WaBeI (Ie,Ba)");
    dpd_buf4_init(&WAbEi, CC3_HC1ET1, C_irr, 27, 29, 27, 29, 0, "Ht_WAbEi (iE,bA)");
    dpd_buf4_init(&WMBIJ, CC3_HC1ET1, C_irr, 0, 20, 2, 20, 0, "Ht_WMBIJ (I>J,MB)");
    dpd_buf4_init(&WMbIj, CC3_HC1ET1, C_irr, 22, 24, 22, 24, 0, "Ht_WMbIj (Ij,Mb)");
    dpd_buf4_init(&WmBiJ, CC3_HC1ET1, C_irr, 23, 27, 23, 27, 0, "Ht_WmBiJ (iJ,mB)");

    dpd_buf4_init(&DIJAB_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    dpd_buf4_init(&DIjAb, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&WAMEF, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    dpd_buf4_init(&WaMeF, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");
    dpd_buf4_init(&WAmEf, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");
    dpd_buf4_init(&WMNIE, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");
    dpd_buf4_init(&WMnIe, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_init(&WmNiE, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");

         /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */

    cc3_sigma_UHF_AAB(&TIJAB, &TIjAb, &TiJaB, &WABEI, &WaBeI, &WAbEi,
       &WMBIJ, &WMbIj, &WmBiJ, 1,  &DIJAB_anti, &DIjAb, &SIA, &Sia,
       1, &FME, &Fme, &WAMEF, &WaMeF, &WAmEf, &WMNIE, &WMnIe, &WmNiE,
       &SIJAB, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
       moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
       moinfo.bvir_off, omega, outfile);

    dpd_buf4_close(&TIJAB); dpd_buf4_close(&TIjAb); dpd_buf4_close(&TiJaB);
    dpd_buf4_close(&WABEI); dpd_buf4_close(&WaBeI); dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&WMBIJ); dpd_buf4_close(&WMbIj); dpd_buf4_close(&WmBiJ);
    dpd_buf4_close(&DIJAB_anti); dpd_buf4_close(&DIjAb);
    dpd_file2_close(&FME); dpd_file2_close(&Fme);
    dpd_buf4_close(&WAMEF); dpd_buf4_close(&WaMeF); dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&WMNIE); dpd_buf4_close(&WMnIe); dpd_buf4_close(&WmNiE);

    /* beta-beta-alpha term 2 */ 

    dpd_buf4_init(&Tijab, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&TIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&TiJaB, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

    dpd_buf4_init(&Wabei, CC3_HC1ET1, C_irr, 30, 15, 30, 17, 0, "Ht_Wabei (ie,b>a)");
    dpd_buf4_init(&WaBeI, CC3_HC1ET1, C_irr, 24, 28, 24, 28, 0, "Ht_WaBeI (Ie,Ba)");
    dpd_buf4_init(&WAbEi, CC3_HC1ET1, C_irr, 27, 29, 27, 29, 0, "Ht_WAbEi (iE,bA)");
    dpd_buf4_init(&Wmbij, CC3_HC1ET1, C_irr, 10, 30, 12, 30, 0, "Ht_Wmbij (i>j,mb)");
    dpd_buf4_init(&WMbIj, CC3_HC1ET1, C_irr, 22, 24, 22, 24, 0, "Ht_WMbIj (Ij,Mb)");
    dpd_buf4_init(&WmBiJ, CC3_HC1ET1, C_irr, 23, 27, 23, 27, 0, "Ht_WmBiJ (iJ,mB)");

    dpd_buf4_init(&Dijab_anti, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    dpd_buf4_init(&DiJaB, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&Wamef, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    dpd_buf4_init(&WaMeF, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");
    dpd_buf4_init(&WAmEf, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");
    dpd_buf4_init(&Wmnie, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");
    dpd_buf4_init(&WMnIe, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_init(&WmNiE, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");

         /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */

    cc3_sigma_UHF_BBA(&Tijab, &TIjAb, &TiJaB, &Wabei, &WaBeI, &WAbEi,
      &Wmbij, &WMbIj, &WmBiJ, 1, &Dijab_anti, &DiJaB, &SIA, &Sia,
      1, &FME, &Fme, &Wamef, &WaMeF, &WAmEf, &Wmnie, &WMnIe, &WmNiE,
      &Sijab, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
      moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
      moinfo.bvir_off, omega, outfile);

    dpd_buf4_close(&Tijab); dpd_buf4_close(&TIjAb); dpd_buf4_close(&TiJaB);
    dpd_buf4_close(&Wabei); dpd_buf4_close(&WaBeI); dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&Wmbij); dpd_buf4_close(&WMbIj); dpd_buf4_close(&WmBiJ);
    dpd_buf4_close(&Dijab_anti); dpd_buf4_close(&DiJaB);
    dpd_file2_close(&FME); dpd_file2_close(&Fme);
    dpd_buf4_close(&Wamef); dpd_buf4_close(&WaMeF); dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&Wmnie); dpd_buf4_close(&WMnIe); dpd_buf4_close(&WmNiE);

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

    dpd_buf4_init(&TIJAB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&WABEI, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    dpd_buf4_init(&WMBIJ, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    dpd_file2_init(&FME, CC3_HC1, C_irr, 0, 1, "HC1 FME");
    dpd_buf4_init(&WAMEF, CC3_HC1, C_irr, 20, 5, 20, 7, 0, "HC1 WAMEF (MA,F>E)");
    dpd_buf4_init(&WMNIE, CC3_HC1, C_irr, 0, 20, 2, 20, 0, "HC1 WMNIE (M>N,IE)");

         /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */

    cc3_sigma_UHF_AAA(&TIJAB, &WABEI, &WMBIJ, 0, NULL, NULL,
        1, &FME, &WAMEF, &WMNIE, &SIJAB, moinfo.aoccpi, moinfo.aocc_off,
        moinfo.avirtpi, moinfo.avir_off, 0.0, outfile);

    dpd_buf4_close(&TIJAB);
    dpd_buf4_close(&WABEI);
    dpd_buf4_close(&WMBIJ);
    dpd_file2_close(&FME);
    dpd_buf4_close(&WAMEF);
    dpd_buf4_close(&WMNIE);

    /*** beta-beta-beta term 3 */

    dpd_buf4_init(&Tijab, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&Wabei, CC3_HET1,0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    dpd_buf4_init(&Wmbij, CC3_HET1,0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    dpd_file2_init(&Fme, CC3_HC1, C_irr, 2, 3, "HC1 Fme");
    dpd_buf4_init(&Wamef, CC3_HC1, C_irr, 30, 15, 30, 17, 0, "HC1 Wamef (ma,f>e)");
    dpd_buf4_init(&Wmnie, CC3_HC1, C_irr, 10, 30, 12, 30, 0, "HC1 Wmnie (m>n,ie)");

         /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */

    cc3_sigma_UHF_BBB(&Tijab, &Wabei, &Wmbij, 0, NULL, NULL,
        1, &Fme, &Wamef, &Wmnie, &Sijab, moinfo.boccpi, moinfo.bocc_off, 
        moinfo.bvirtpi, moinfo.bvir_off, 0.0, outfile);

    dpd_buf4_close(&Tijab);
    dpd_buf4_close(&Wabei);
    dpd_buf4_close(&Wmbij);
    dpd_file2_close(&Fme);
    dpd_buf4_close(&Wamef);
    dpd_buf4_close(&Wmnie);

    /*** alpha-alpha-beta term 3 */

    dpd_buf4_init(&TIJAB, CC_TAMPS, 0,  0,  5,  2,  7, 0, "tIJAB");
    dpd_buf4_init(&TIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&TiJaB, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

    dpd_buf4_init(&WABEI, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    dpd_buf4_init(&WaBeI, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
    dpd_buf4_init(&WAbEi, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
    dpd_buf4_init(&WMBIJ, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    dpd_buf4_init(&WMbIj, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
    dpd_buf4_init(&WmBiJ, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

    dpd_file2_init(&FME, CC3_HC1, C_irr, 0, 1, "HC1 FME");
    dpd_file2_init(&Fme, CC3_HC1, C_irr, 2, 3, "HC1 Fme");
    dpd_buf4_init(&WAMEF, CC3_HC1, C_irr, 20, 5, 20, 7, 0, "HC1 WAMEF (MA,F>E)");
    dpd_buf4_init(&WaMeF, CC3_HC1, C_irr, 24, 28, 24, 28, 0, "HC1 WaMeF (Ma,Fe)");
    dpd_buf4_init(&WAmEf, CC3_HC1, C_irr, 27, 29, 27, 29, 0, "HC1 WAmEf (mA,fE)");
    dpd_buf4_init(&WMNIE, CC3_HC1, C_irr, 0, 20, 2, 20, 0, "HC1 WMNIE (M>N,IE)");
    dpd_buf4_init(&WMnIe, CC3_HC1, C_irr, 22, 24, 22, 24, 0, "HC1 WMnIe (Mn,Ie)");
    dpd_buf4_init(&WmNiE, CC3_HC1, C_irr, 23, 27, 23, 27, 0, "HC1 WmNiE (mN,iE)");

         /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */

    cc3_sigma_UHF_AAB(&TIJAB, &TIjAb, &TiJaB, &WABEI, &WaBeI, &WAbEi,
       &WMBIJ, &WMbIj, &WmBiJ, 0, NULL, NULL, NULL, NULL,
       1, &FME, &Fme, &WAMEF, &WaMeF, &WAmEf, &WMNIE, &WMnIe, &WmNiE,
       &SIJAB, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
       moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
       moinfo.bvir_off, 0.0, outfile);

    dpd_buf4_close(&TIJAB); dpd_buf4_close(&TIjAb); dpd_buf4_close(&TiJaB);
    dpd_buf4_close(&WABEI); dpd_buf4_close(&WaBeI); dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&WMBIJ); dpd_buf4_close(&WMbIj); dpd_buf4_close(&WmBiJ);
    dpd_file2_close(&FME); dpd_file2_close(&Fme);
    dpd_buf4_close(&WAMEF); dpd_buf4_close(&WaMeF); dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&WMNIE); dpd_buf4_close(&WMnIe); dpd_buf4_close(&WmNiE);

    /*** beta-beta-alpha term 3 */ 

    dpd_buf4_init(&Tijab, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&TIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&TiJaB, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

    dpd_buf4_init(&Wabei, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    dpd_buf4_init(&WaBeI, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
    dpd_buf4_init(&WAbEi, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
    dpd_buf4_init(&Wmbij, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    dpd_buf4_init(&WMbIj, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
    dpd_buf4_init(&WmBiJ, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

    dpd_file2_init(&FME, CC3_HC1, C_irr, 0, 1, "HC1 FME");
    dpd_file2_init(&Fme, CC3_HC1, C_irr, 2, 3, "HC1 Fme");
    dpd_buf4_init(&Wamef, CC3_HC1, C_irr, 30, 15, 30, 17, 0, "HC1 Wamef (ma,f>e)");
    dpd_buf4_init(&WaMeF, CC3_HC1, C_irr, 24, 28, 24, 28, 0, "HC1 WaMeF (Ma,Fe)");
    dpd_buf4_init(&WAmEf, CC3_HC1, C_irr, 27, 29, 27, 29, 0, "HC1 WAmEf (mA,fE)");
    dpd_buf4_init(&Wmnie, CC3_HC1, C_irr, 10, 30, 12, 30, 0, "HC1 Wmnie (m>n,ie)");
    dpd_buf4_init(&WMnIe, CC3_HC1, C_irr, 22, 24, 22, 24, 0, "HC1 WMnIe (Mn,Ie)");
    dpd_buf4_init(&WmNiE, CC3_HC1, C_irr, 23, 27, 23, 27, 0, "HC1 WmNiE (mN,iE)");

         /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */

    cc3_sigma_UHF_BBA(&Tijab, &TIjAb, &TiJaB, &Wabei, &WaBeI, &WAbEi,
      &Wmbij, &WMbIj, &WmBiJ, 0, NULL, NULL, NULL, NULL,
      1, &FME, &Fme, &Wamef, &WaMeF, &WAmEf, &Wmnie, &WMnIe, &WmNiE,
      &Sijab, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
      moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
      moinfo.bvir_off, 0.0, outfile);

    dpd_buf4_close(&Tijab); dpd_buf4_close(&TIjAb); dpd_buf4_close(&TiJaB);
    dpd_buf4_close(&Wabei); dpd_buf4_close(&WaBeI); dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&Wmbij); dpd_buf4_close(&WMbIj); dpd_buf4_close(&WmBiJ);
    dpd_file2_close(&FME); dpd_file2_close(&Fme);
    dpd_buf4_close(&Wamef); dpd_buf4_close(&WaMeF); dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&Wmnie); dpd_buf4_close(&WMnIe); dpd_buf4_close(&WmNiE);

    dpd_file2_close(&SIA);
    dpd_file2_close(&Sia);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
  check_sum("<Psi|H'<T|(Uhat T2)c|0>|T>/(w-wt)", i, C_irr);
#endif
  }
  return;
}

}} // namespace psi::cceom
