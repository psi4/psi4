/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/* cc3(): calls libdpd functions to compute the triples contributions
 * to t1 and t2,
 *  <S| H    <T| (Uhat T2)c   |0> |T> -> t1
 *  <D| Hhat <T| (Uhat T2)c   |0> |T> -> t2 
 * Denominators are subsequently applied in denom()
*/

void cc3(void)
{
  dpdfile2 TIA_new, Tia_new, FME, Fme;
  dpdbuf4 TIJAB_new, Tijab_new, TIjAb_new;
  dpdbuf4 TIJAB, Tijab, TIjAb, TiJaB;
  dpdbuf4 WABEI, WAbEi, Wabei, WaBeI;
  dpdbuf4 WMBIJ, WMbIj, Wmbij, WmBiJ;
  dpdbuf4 Dints, DIJAB_anti, Dijab_anti, DIjAb, DiJaB;
  dpdbuf4 WAMEF, WAmEf, Wamef, WaMeF;
  dpdbuf4 WMNIE, WMnIe, Wmnie, WmNiE;

  if(params.ref == 0) { /* RHF */
    dpd_file2_init(&TIA_new, CC_OEI, 0, 0, 1, "New tIA");
    dpd_buf4_init(&TIjAb_new, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_buf4_init(&TIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&WAbEi, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
    dpd_buf4_init(&WMbIj, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
    dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&WAmEf, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAmEf (mA,Ef)");
    dpd_buf4_init(&WMnIe, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");

    if (params.t3_Ws_incore)
      cc3_sigma_RHF_ic(&TIjAb, &WAbEi, &WMbIj, 1, &Dints, &TIA_new, 1, &FME, &WAmEf,
        &WMnIe, &TIjAb_new, moinfo.occpi, moinfo.occ_off, moinfo.virtpi,
        moinfo.vir_off, 0.0, outfile, params.nthreads, params.newtrips);
    else 
      cc3_sigma_RHF(&TIjAb, &WAbEi, &WMbIj, 1, &Dints, &TIA_new, 1, &FME, &WAmEf,
        &WMnIe, &TIjAb_new, moinfo.occpi, moinfo.occ_off, moinfo.virtpi,
        moinfo.vir_off, 0.0, outfile, params.newtrips);

    dpd_buf4_close(&TIjAb);
    dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&WMbIj);
    dpd_buf4_close(&Dints);
    dpd_file2_close(&FME);
    dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&WMnIe);

    dpd_file2_close(&TIA_new);
    dpd_buf4_close(&TIjAb_new);
  }
  else if(params.ref == 2) { /* UHF */
    dpd_file2_init(&TIA_new, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&Tia_new, CC_OEI, 0, 2, 3, "New tia");
    dpd_buf4_init(&TIJAB_new, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&Tijab_new, CC_TAMPS, 0, 10, 15, 12, 17, 0, "New tijab");
    dpd_buf4_init(&TIjAb_new, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    /*** alpha-alpha-alpha */

    dpd_buf4_init(&TIJAB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&WABEI, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    dpd_buf4_init(&WMBIJ, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    dpd_buf4_init(&DIJAB_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&WAMEF, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    dpd_buf4_init(&WMNIE, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");

    cc3_sigma_UHF_AAA(&TIJAB, &WABEI, &WMBIJ, 1, &DIJAB_anti, &TIA_new,
        1, &FME, &WAMEF, &WMNIE, &TIJAB_new, moinfo.aoccpi, moinfo.aocc_off,
        moinfo.avirtpi, moinfo.avir_off, 0.0, outfile);

    dpd_buf4_close(&TIJAB);
    dpd_buf4_close(&WABEI);
    dpd_buf4_close(&WMBIJ);
    dpd_buf4_close(&DIJAB_anti);
    dpd_file2_close(&FME);
    dpd_buf4_close(&WAMEF);
    dpd_buf4_close(&WMNIE);

    /*** beta-beta-beta */

    dpd_buf4_init(&Tijab, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&Wabei, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    dpd_buf4_init(&Wmbij, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    dpd_buf4_init(&Dijab_anti, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&Wamef, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    dpd_buf4_init(&Wmnie, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");

    cc3_sigma_UHF_BBB(&Tijab, &Wabei, &Wmbij, 1, &Dijab_anti, &Tia_new,
        1, &Fme, &Wamef, &Wmnie, &Tijab_new, moinfo.boccpi, moinfo.bocc_off,
        moinfo.bvirtpi, moinfo.bvir_off, 0.0, outfile);

    dpd_buf4_close(&Tijab);
    dpd_buf4_close(&Wabei);
    dpd_buf4_close(&Wmbij);
    dpd_buf4_close(&Dijab_anti);
    dpd_file2_close(&Fme);
    dpd_buf4_close(&Wamef);
    dpd_buf4_close(&Wmnie);

    /*** alpha-alpha-beta */ 

    dpd_buf4_init(&TIJAB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&TIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&TiJaB, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

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

    cc3_sigma_UHF_AAB(&TIJAB, &TIjAb, &TiJaB, &WABEI, &WaBeI, &WAbEi,
       &WMBIJ, &WMbIj, &WmBiJ, 1,  &DIJAB_anti, &DIjAb, &TIA_new, &Tia_new,
       1, &FME, &Fme, &WAMEF, &WaMeF, &WAmEf, &WMNIE, &WMnIe, &WmNiE,
       &TIJAB_new, &TIjAb_new, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
       moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
       moinfo.bvir_off, 0.0, outfile);

    dpd_buf4_close(&TIJAB); dpd_buf4_close(&TIjAb); dpd_buf4_close(&TiJaB);
    dpd_buf4_close(&WABEI); dpd_buf4_close(&WaBeI); dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&WMBIJ); dpd_buf4_close(&WMbIj); dpd_buf4_close(&WmBiJ);
    dpd_buf4_close(&DIJAB_anti); dpd_buf4_close(&DIjAb);
    dpd_file2_close(&FME); dpd_file2_close(&Fme);
    dpd_buf4_close(&WAMEF); dpd_buf4_close(&WaMeF); dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&WMNIE); dpd_buf4_close(&WMnIe); dpd_buf4_close(&WmNiE);

    /*** beta-beta-alpha term 1 */ 

    dpd_buf4_init(&Tijab, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&TIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&TiJaB, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

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

    cc3_sigma_UHF_BBA(&Tijab, &TIjAb, &TiJaB, &Wabei, &WaBeI, &WAbEi,
      &Wmbij, &WMbIj, &WmBiJ, 1, &Dijab_anti, &DiJaB, &TIA_new, &Tia_new,
      1, &FME, &Fme, &Wamef, &WaMeF, &WAmEf, &Wmnie, &WMnIe, &WmNiE,
      &Tijab_new, &TIjAb_new, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
      moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
      moinfo.bvir_off, 0.0, outfile);

    dpd_buf4_close(&Tijab); dpd_buf4_close(&TIjAb); dpd_buf4_close(&TiJaB);
    dpd_buf4_close(&Wabei); dpd_buf4_close(&WaBeI); dpd_buf4_close(&WAbEi);
    dpd_buf4_close(&Wmbij); dpd_buf4_close(&WMbIj); dpd_buf4_close(&WmBiJ);
    dpd_buf4_close(&Dijab_anti); dpd_buf4_close(&DiJaB);
    dpd_file2_close(&FME); dpd_file2_close(&Fme);
    dpd_buf4_close(&Wamef); dpd_buf4_close(&WaMeF); dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&Wmnie); dpd_buf4_close(&WMnIe); dpd_buf4_close(&WmNiE);

    /*
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_print(&t1, outfile);
    dpd_file2_close(&t1);
    dpd_file2_init(&t1, CC_OEI, 0, 2, 3, "New tia");
    dpd_file2_print(&t1, outfile);
    dpd_file2_close(&t1);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_print(&T2, outfile, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_buf4_print(&T2, outfile, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_buf4_print(&T2, outfile, 1);
    dpd_buf4_close(&T2); 
    */
  }
}
}} // namespace psi::ccenergy
