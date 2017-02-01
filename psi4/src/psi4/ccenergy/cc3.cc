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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "Params.h"
#include "MOInfo.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/* cc3(): calls libdpd functions to compute the triples contributions
 * to t1 and t2,
 *  <S| H    <T| (Uhat T2)c   |0> |T> -> t1
 *  <D| Hhat <T| (Uhat T2)c   |0> |T> -> t2
 * Denominators are subsequently applied in denom()
*/

void CCEnergyWavefunction::cc3(void)
{
  dpdfile2 TIA_new, Tia_new, FME, Fme;
  dpdbuf4 TIJAB_new, Tijab_new, TIjAb_new;
  dpdbuf4 TIJAB, Tijab, TIjAb, TiJaB;
  dpdbuf4 WABEI, WAbEi, Wabei, WaBeI;
  dpdbuf4 WMBIJ, WMbIj, Wmbij, WmBiJ;
  dpdbuf4 Dints, DIJAB_anti, Dijab_anti, DIjAb, DiJaB;
  dpdbuf4 WAMEF, WAmEf, Wamef, WaMeF;
  dpdbuf4 WMNIE, WMnIe, Wmnie, WmNiE;

  if(params_.ref == 0) { /* RHF */
    global_dpd_->file2_init(&TIA_new, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->buf4_init(&TIjAb_new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
    global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAmEf (mA,Ef)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");

    if (params_.t3_Ws_incore)
      global_dpd_->cc3_sigma_RHF_ic(&TIjAb, &WAbEi, &WMbIj, 1, &Dints, &TIA_new, 1, &FME, &WAmEf,
        &WMnIe, &TIjAb_new, moinfo_.occpi, moinfo_.occ_off, moinfo_.virtpi,
        moinfo_.vir_off, 0.0, "outfile", params_.nthreads, params_.newtrips);
    else
      global_dpd_->cc3_sigma_RHF(&TIjAb, &WAbEi, &WMbIj, 1, &Dints, &TIA_new, 1, &FME, &WAmEf,
        &WMnIe, &TIjAb_new, moinfo_.occpi, moinfo_.occ_off, moinfo_.virtpi,
        moinfo_.vir_off, 0.0, "outfile", params_.newtrips);

    global_dpd_->buf4_close(&TIjAb);
    global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&WMbIj);
    global_dpd_->buf4_close(&Dints);
    global_dpd_->file2_close(&FME);
    global_dpd_->buf4_close(&WAmEf);
    global_dpd_->buf4_close(&WMnIe);

    global_dpd_->file2_close(&TIA_new);
    global_dpd_->buf4_close(&TIjAb_new);
  }
  else if(params_.ref == 2) { /* UHF */
    global_dpd_->file2_init(&TIA_new, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_init(&Tia_new, PSIF_CC_OEI, 0, 2, 3, "New tia");
    global_dpd_->buf4_init(&TIJAB_new, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&Tijab_new, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "New tijab");
    global_dpd_->buf4_init(&TIjAb_new, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    /*** alpha-alpha-alpha */

    global_dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&WABEI, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    global_dpd_->buf4_init(&DIJAB_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");

    global_dpd_->cc3_sigma_UHF_AAA(&TIJAB, &WABEI, &WMBIJ, 1, &DIJAB_anti, &TIA_new,
        1, &FME, &WAMEF, &WMNIE, &TIJAB_new, moinfo_.aoccpi, moinfo_.aocc_off,
        moinfo_.avirtpi, moinfo_.avir_off, 0.0, "outfile");

    global_dpd_->buf4_close(&TIJAB);
    global_dpd_->buf4_close(&WABEI);
    global_dpd_->buf4_close(&WMBIJ);
    global_dpd_->buf4_close(&DIJAB_anti);
    global_dpd_->file2_close(&FME);
    global_dpd_->buf4_close(&WAMEF);
    global_dpd_->buf4_close(&WMNIE);

    /*** beta-beta-beta */

    global_dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&Wabei, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    global_dpd_->buf4_init(&Wmbij, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    global_dpd_->buf4_init(&Dijab_anti, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");

    global_dpd_->cc3_sigma_UHF_BBB(&Tijab, &Wabei, &Wmbij, 1, &Dijab_anti, &Tia_new,
        1, &Fme, &Wamef, &Wmnie, &Tijab_new, moinfo_.boccpi, moinfo_.bocc_off,
        moinfo_.bvirtpi, moinfo_.bvir_off, 0.0, "outfile");

    global_dpd_->buf4_close(&Tijab);
    global_dpd_->buf4_close(&Wabei);
    global_dpd_->buf4_close(&Wmbij);
    global_dpd_->buf4_close(&Dijab_anti);
    global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&Wamef);
    global_dpd_->buf4_close(&Wmnie);

    /*** alpha-alpha-beta */

    global_dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&TiJaB, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

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

    global_dpd_->cc3_sigma_UHF_AAB(&TIJAB, &TIjAb, &TiJaB, &WABEI, &WaBeI, &WAbEi,
       &WMBIJ, &WMbIj, &WmBiJ, 1,  &DIJAB_anti, &DIjAb, &TIA_new, &Tia_new,
       1, &FME, &Fme, &WAMEF, &WaMeF, &WAmEf, &WMNIE, &WMnIe, &WmNiE,
       &TIJAB_new, &TIjAb_new, moinfo_.aoccpi, moinfo_.aocc_off, moinfo_.boccpi,
       moinfo_.bocc_off, moinfo_.avirtpi, moinfo_.avir_off, moinfo_.bvirtpi,
       moinfo_.bvir_off, 0.0, "outfile");

    global_dpd_->buf4_close(&TIJAB); global_dpd_->buf4_close(&TIjAb); global_dpd_->buf4_close(&TiJaB);
    global_dpd_->buf4_close(&WABEI); global_dpd_->buf4_close(&WaBeI); global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&WMBIJ); global_dpd_->buf4_close(&WMbIj); global_dpd_->buf4_close(&WmBiJ);
    global_dpd_->buf4_close(&DIJAB_anti); global_dpd_->buf4_close(&DIjAb);
    global_dpd_->file2_close(&FME); global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&WAMEF); global_dpd_->buf4_close(&WaMeF); global_dpd_->buf4_close(&WAmEf);
    global_dpd_->buf4_close(&WMNIE); global_dpd_->buf4_close(&WMnIe); global_dpd_->buf4_close(&WmNiE);

    /*** beta-beta-alpha term 1 */

    global_dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&TiJaB, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

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

    global_dpd_->cc3_sigma_UHF_BBA(&Tijab, &TIjAb, &TiJaB, &Wabei, &WaBeI, &WAbEi,
      &Wmbij, &WMbIj, &WmBiJ, 1, &Dijab_anti, &DiJaB, &TIA_new, &Tia_new,
      1, &FME, &Fme, &Wamef, &WaMeF, &WAmEf, &Wmnie, &WMnIe, &WmNiE,
      &Tijab_new, &TIjAb_new, moinfo_.aoccpi, moinfo_.aocc_off, moinfo_.boccpi,
      moinfo_.bocc_off, moinfo_.avirtpi, moinfo_.avir_off, moinfo_.bvirtpi,
      moinfo_.bvir_off, 0.0, "outfile");

    global_dpd_->buf4_close(&Tijab); global_dpd_->buf4_close(&TIjAb); global_dpd_->buf4_close(&TiJaB);
    global_dpd_->buf4_close(&Wabei); global_dpd_->buf4_close(&WaBeI); global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&Wmbij); global_dpd_->buf4_close(&WMbIj); global_dpd_->buf4_close(&WmBiJ);
    global_dpd_->buf4_close(&Dijab_anti); global_dpd_->buf4_close(&DiJaB);
    global_dpd_->file2_close(&FME); global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&Wamef); global_dpd_->buf4_close(&WaMeF); global_dpd_->buf4_close(&WAmEf);
    global_dpd_->buf4_close(&Wmnie); global_dpd_->buf4_close(&WMnIe); global_dpd_->buf4_close(&WmNiE);

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
