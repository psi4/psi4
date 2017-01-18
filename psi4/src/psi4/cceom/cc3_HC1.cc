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
#include <cstdlib>
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void HC1_F(int i, int C_irr);
void HC1_Wamef(int i, int C_irr);
void HC1_Wmnie(int i, int C_irr);
void HC1_Wmnij(int i, int C_irr);
void HC1_Wmbij(int i, int C_irr);
void HC1_Wmbej(int i, int C_irr);
void HC1_Wabei(int i, int C_irr);
void purge_HC1(int C_irr);
void purge_Wmnie(int C_irr);
void purge_Wmbij(int C_irr);
void purge_Wabei(int C_irr);

void cc3_HC1 (int i, int C_irr) {
  HC1_F(i, C_irr);
  HC1_Wamef(i, C_irr);
  HC1_Wmnie(i, C_irr);
/*  HC1_Wmnij(i, C_irr); */
/*  HC1_Wmbij(i, C_irr); */
/*  HC1_Wmbej(i, C_irr); */
  HC1_Wabei(i, C_irr);

  if (params.ref == 1) purge_HC1(C_irr);

  return;
}

/* constructs matrix elements of [H, C1] for CC3 EOM code */

void HC1_F(int i, int C_irr) {

  dpdfile2 FME, Fme, Cme, CME;
  dpdbuf4 D_anti, D;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);
	double tval;

  if(params.ref == 0) {
    /* HC1_F()  Fme = +C_n^f <mn||ef> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);

    global_dpd_->file2_init(&FME, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 FME");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->dot13(&CME, &D, &FME, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&FME);

    global_dpd_->file2_close(&CME);
  }

  else if(params.ref == 1) { /** ROHF **/
    /* HC1_F()  Fme = +C_n^f <mn||ef> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);

    global_dpd_->file2_init(&FME, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 FME");
    global_dpd_->file2_init(&Fme, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 Fme");

    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->dot13(&CME, &D_anti, &FME, 0, 0, 1.0, 0.0);
    global_dpd_->dot13(&Cme, &D,      &FME, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&Cme, &D_anti, &Fme, 0, 0, 1.0, 0.0);
    global_dpd_->dot13(&CME, &D,      &Fme, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D_anti);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&Fme);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/
    /* HC1_F()  Fme = +C_n^f <mn||ef> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);

    global_dpd_->file2_init(&FME, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 FME");
    global_dpd_->file2_init(&Fme, PSIF_CC3_HC1, C_irr, 2, 3, "HC1 Fme");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->contract422(&D, &CME, &FME, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->contract422(&D, &Cme, &FME, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->contract422(&D, &Cme, &Fme, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->contract422(&D, &CME, &Fme, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&Fme);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }

  return;
}


void HC1_Wamef(int i, int C_irr) {
  dpdbuf4 Wamef, WAMEF, WAmEf, WaMeF, W, D_a, D;
  dpdfile2 Cme, CME;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if(params.ref == 0) { /** RHF **/
    /* HC1_Wamef():  Wamef = -Cna <nm||ef> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 11, 5, 11, 5, 0, "HC1 WAmEf (Am,Ef)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&CME, &D, &W, 0, 0, 0, -1, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, qprs, 10, 5, "HC1 WAmEf (mA,Ef)");
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&CME);
  }

  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);

    /* C(N,A) <NM||EF> --> W(AM,E>F) */
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HC1, C_irr, 11, 7, 11, 7, 0, "HC1 WAMEF (AM,E>F)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&CME, &D, &WAMEF, 0, 0, 0, -1, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WAMEF);

    /* C(n,a) <nm||ef> --> W(am,e>f) */
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HC1, C_irr, 11, 7, 11, 7, 0, "HC1 Wamef (am,e>f)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&Cme, &D, &Wamef, 0, 0, 0, -1, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Wamef);

    /* C(N,A) <Nm|Ef> --> W(Am,Ef) */
    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HC1, C_irr, 11, 5, 11, 5, 0, "HC1 WAmEf (Am,Ef)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&CME, &D, &WAmEf, 0, 0, 0, -1, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WAmEf);

    /* C(n,a) <nM|eF> --> W(aM,eF) */
    global_dpd_->buf4_init(&WaMeF, PSIF_CC3_HC1, C_irr, 11, 5, 11, 5, 0, "HC1 WaMeF (aM,eF)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&Cme, &D, &WaMeF, 0, 0, 0, -1, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WaMeF);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);

  } /** ROHF **/
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /* T(N,A) <NM||EF> --> W(AM,E>F) */
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HC1, C_irr, 21, 7, 21, 7, 0, "HC1 WAMEF (AM,E>F)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    global_dpd_->contract244(&CME, &D, &WAMEF, 0, 0, 0, -1, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WAMEF);

    /* T(n,a) <nm||ef> --> W(am,e>f) */
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HC1, C_irr, 31, 17, 31, 17, 0, "HC1 Wamef (am,e>f)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&Cme, &D, &Wamef, 0, 0, 0, -1, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Wamef);

    /* T(N,A) <Nm|Ef> --> W(Am,Ef) */
    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HC1, C_irr, 26, 28, 26, 28, 0, "HC1 WAmEf (Am,Ef)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract244(&CME, &D, &WAmEf, 0, 0, 0, -1, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WAmEf);

    /* T(n,a) <nM|eF> --> W(aM,eF) */
    global_dpd_->buf4_init(&WaMeF, PSIF_CC3_HC1, C_irr, 25, 29, 25, 29, 0, "HC1 WaMeF (aM,eF)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract244(&Cme, &D, &WaMeF, 0, 0, 0, -1, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WaMeF);

    /* Do final sort */
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HC1, C_irr, 21, 7, 21, 7, 0, "HC1 WAMEF (AM,E>F)");
    global_dpd_->buf4_sort(&WAMEF, PSIF_CC3_HC1, qprs, 20, 7, "HC1 WAMEF (MA,F>E)");
    global_dpd_->buf4_close(&WAMEF);
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HC1, C_irr, 20, 7, 20, 7, 0, "HC1 WAMEF (MA,F>E)");
    global_dpd_->buf4_scm(&WAMEF, -1.0);
    global_dpd_->buf4_close(&WAMEF);

    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HC1, C_irr, 31, 17, 31, 17, 0, "HC1 Wamef (am,e>f)");
    global_dpd_->buf4_sort(&Wamef, PSIF_CC3_HC1, qprs, 30, 17, "HC1 Wamef (ma,f>e)");
    global_dpd_->buf4_close(&Wamef);
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HC1, C_irr, 30, 17, 30, 17, 0, "HC1 Wamef (ma,f>e)");
    global_dpd_->buf4_scm(&Wamef, -1.0);
    global_dpd_->buf4_close(&Wamef);

    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HC1, C_irr, 26, 28, 26, 28, 0, "HC1 WAmEf (Am,Ef)");
    global_dpd_->buf4_sort(&WAmEf, PSIF_CC3_HC1, qpsr, 27, 29, "HC1 WAmEf (mA,fE)");
    global_dpd_->buf4_close(&WAmEf);

    global_dpd_->buf4_init(&WaMeF, PSIF_CC3_HC1, C_irr, 25, 29, 25, 29, 0, "HC1 WaMeF (aM,eF)");
    global_dpd_->buf4_sort(&WaMeF, PSIF_CC3_HC1, qpsr, 24, 28, "HC1 WaMeF (Ma,Fe)");
    global_dpd_->buf4_close(&WaMeF);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  } /** UHF **/

  return;
}


void HC1_Wmnie(int i, int C_irr) {
  dpdbuf4 W, Wmnie, WMNIE, WMnIe, WmNiE, WMniE, WmNIe;
  dpdbuf4 D, D_a, Z;
  dpdfile2 Cme, CME;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if(params.ref == 0) { /** RHF **/
    /* HC1_Wmnie():  Wmnie = + Cif <mn||fe> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);

    /* C(I,F) * D(Mn,Fe) --> W(Mn,Ie) */
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HC1, C_irr, 0, 10, 0, 10, 0, "HC1 WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&CME, &D, &WMnIe, 1, 2, 1, 1, 0.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&D);
    /* W(Mn,Ie) --> W(Mn,eI) */
    /* dpd_buf4_sort(&WMnIe, CC3_HC1, pqsr, 0, 11, "HC1 WMnIe (Mn,eI)"); */
    global_dpd_->buf4_close(&WMnIe);
  }

  else if(params.ref == 1) { /** ROHF **/
    /* HC1_Wmnie():  Wmnie = + Cif <mn||fe> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);

    /* D(M>N,EF) * C(I,F) --> W(M>N,EI) */
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HC1, C_irr, 2, 11, 2, 11, 0, "HC1 WMNIE (M>N,EI)");
    global_dpd_->buf4_init(&D_a, PSIF_CC_DINTS, 0, 2, 5, 2, 5,0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D_a,&CME,&WMNIE, 3, 1, 0, -1, 0);
    global_dpd_->buf4_close(&D_a);
    global_dpd_->buf4_close(&WMNIE);

    /* D(m>n,ef) * C(i,f) --> W(m>n,ei) */
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HC1, C_irr, 2, 11, 2, 11, 0, "HC1 Wmnie (m>n,ei)");
    global_dpd_->buf4_init(&D_a, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D_a, &Cme, &Wmnie, 3, 1, 0, -1, 0);
    global_dpd_->buf4_close(&D_a);
    global_dpd_->buf4_close(&Wmnie);

    /* D(Mn,Fe) * C(I,F) --> W(Mn,Ie) */
    global_dpd_->buf4_init(&WMnIe, PSIF_CC_TMP0, C_irr, 0, 10, 0, 10, 0, "HC1 WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&CME, &D, &WMnIe, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    /* W(Mn,Ie) --> W(Mn,eI) */
    global_dpd_->buf4_sort(&WMnIe, PSIF_CC3_HC1, pqsr, 0, 11, "HC1 WMnIe (Mn,eI)");
    global_dpd_->buf4_close(&WMnIe);

    /* D(mN,fE) * C(i,f) --> W(mN,iE) */
    global_dpd_->buf4_init(&WmNiE, PSIF_CC_TMP1, C_irr, 0, 10, 0, 10, 0, "HC1 WmNiE (mN,iE)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&Cme,&D,&WmNiE, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    /* W(mN,iE) --> W(mN,Ei) */
    global_dpd_->buf4_sort(&WmNiE, PSIF_CC3_HC1, pqsr, 0, 11, "HC1 WmNiE (mN,Ei)");
    global_dpd_->buf4_close(&WmNiE);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);

    purge_Wmnie(C_irr); /* before sorting here */

    /* also put "normal" sorted versions in CC_HBAR */
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HC1, C_irr, 2, 11, 2, 11, 0, "HC1 WMNIE (M>N,EI)");
    global_dpd_->buf4_sort(&WMNIE, PSIF_CC3_HC1, pqsr, 2, 10, "HC1 WMNIE (M>N,IE)");
    global_dpd_->buf4_close(&WMNIE);
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HC1, C_irr, 2, 11, 2, 11, 0, "HC1 Wmnie (m>n,ei)");
    global_dpd_->buf4_sort(&Wmnie, PSIF_CC3_HC1, pqsr, 2, 10, "HC1 Wmnie (m>n,ie)");
    global_dpd_->buf4_close(&Wmnie);
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HC1, C_irr, 0, 11, 0, 11, 0, "HC1 WMnIe (Mn,eI)");
    global_dpd_->buf4_sort(&WMnIe, PSIF_CC3_HC1, pqsr, 0, 10, "HC1 WMnIe (Mn,Ie)");
    global_dpd_->buf4_close(&WMnIe);
    global_dpd_->buf4_init(&WmNiE, PSIF_CC3_HC1, C_irr, 0, 11, 0, 11, 0, "HC1 WmNiE (mN,Ei)");
    global_dpd_->buf4_sort(&WmNiE, PSIF_CC3_HC1, pqsr, 0, 10, "HC1 WmNiE (mN,iE)");
    global_dpd_->buf4_close(&WmNiE);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /* <M>N||EF> T(I,F) --> W(M>N,EI) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 2, 21, 2, 21, 0, "HC1 WMNIE (M>N,EI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    global_dpd_->contract424(&D, &CME, &W, 3, 1, 0, -1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /* <m>n||ef> T(i,f) --> W(m>n,ei) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 12, 31, 12, 31, 0, "HC1 Wmnie (m>n,ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D, &Cme, &W, 3, 1, 0, -1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /* Z(nM,eI) = <nM|eF> T(I,F) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, C_irr, 23, 25, 23, 25, 0, "Z(nM,eI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract424(&D, &CME, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* Z(nM,eI) --> W(Mn,eI) */
    global_dpd_->buf4_sort(&Z, PSIF_CC3_HC1, qprs, 22, 25, "HC1 WMnIe (Mn,eI)");
    global_dpd_->buf4_close(&Z);

    /* Z(Nm,Ei) = <Nm|Ef> T(i,f) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, C_irr, 22, 26, 22, 26, 0, "Z(Nm,Ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract424(&D, &Cme, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* Z(Nm,Ei) --> W(mN,Ei) */
    global_dpd_->buf4_sort(&Z, PSIF_CC3_HC1, qprs, 23, 26, "HC1 WmNiE (mN,Ei)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);

    /* also put "normal" sorted versions in CC3_HC1 */
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HC1, C_irr, 2, 21, 2, 21, 0, "HC1 WMNIE (M>N,EI)");
    global_dpd_->buf4_sort(&WMNIE, PSIF_CC3_HC1, pqsr, 2, 20, "HC1 WMNIE (M>N,IE)");
    global_dpd_->buf4_close(&WMNIE);
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HC1, C_irr, 12, 31, 12, 31, 0, "HC1 Wmnie (m>n,ei)");
    global_dpd_->buf4_sort(&Wmnie, PSIF_CC3_HC1, pqsr, 12, 30, "HC1 Wmnie (m>n,ie)");
    global_dpd_->buf4_close(&Wmnie);
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HC1, C_irr, 22, 25, 22, 25, 0, "HC1 WMnIe (Mn,eI)");
    global_dpd_->buf4_sort(&WMnIe, PSIF_CC3_HC1, pqsr, 22, 24, "HC1 WMnIe (Mn,Ie)");
    global_dpd_->buf4_close(&WMnIe);
    global_dpd_->buf4_init(&WmNiE, PSIF_CC3_HC1, C_irr, 23, 26, 23, 26, 0, "HC1 WmNiE (mN,Ei)");
    global_dpd_->buf4_sort(&WmNiE, PSIF_CC3_HC1, pqsr, 23, 27, "HC1 WmNiE (mN,iE)");
    global_dpd_->buf4_close(&WmNiE);
  }

  return;
}

void HC1_Wmnij(int i, int C_irr)
{
  dpdbuf4 WMNIJ, Wmnij, WMnIj, W;
  dpdbuf4 Eijka, Eijka_anti, Eaijk, Eaijk_anti;
  dpdfile2 Cme, CME;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if(params.ref == 0) { /** RHF **/
    /** HC1_Wmnij():  Wmnij = + P(ij) Cje <mn||ie> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);

    global_dpd_->buf4_init(&WMnIj, PSIF_CC3_HC1, C_irr, 0, 0, 0, 0, 0, "HC1 WMnIj (Mn,Ij)");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->contract244(&CME, &Eaijk, &WMnIj, 1, 0, 1, 1, 0.0);
    global_dpd_->buf4_close(&Eaijk);

    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->contract424(&Eijka, &CME, &WMnIj, 3, 1, 0, 1, 1.0);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&WMnIj);

    global_dpd_->file2_close(&CME);
  }

  else if(params.ref == 1) { /** ROHF **/
    /** HC1_Wmnij():  Wmnij = + P(ij) Cje <mn||ie> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);

    global_dpd_->buf4_init(&Eijka_anti, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&Eaijk_anti, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

    global_dpd_->buf4_init(&WMNIJ, PSIF_CC3_HC1, C_irr, 2, 0, 2, 2, 0, "HC1 WMNIJ (M>N,I>J)");
    global_dpd_->contract424(&Eijka_anti, &CME, &WMNIJ, 3, 1, 0, 1, 0.0);
    global_dpd_->contract244(&CME, &Eaijk_anti, &WMNIJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMNIJ);

    global_dpd_->buf4_init(&Wmnij, PSIF_CC3_HC1, C_irr, 2, 0, 2, 2, 0, "HC1 Wmnij (m>n,i>j)");
    global_dpd_->contract424(&Eijka_anti, &Cme, &Wmnij, 3, 1, 0, 1, 0.0);
    global_dpd_->contract244(&Cme, &Eaijk_anti, &Wmnij, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Wmnij);

    global_dpd_->buf4_init(&WMnIj, PSIF_CC3_HC1, C_irr, 0, 0, 0, 0, 0, "HC1 WMnIj (Mn,Ij)");
    global_dpd_->contract424(&Eijka, &Cme, &WMnIj, 3, 1, 0, 1, 0.0);
    global_dpd_->contract244(&CME, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMnIj);

    global_dpd_->buf4_close(&Eijka_anti);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk_anti);
    global_dpd_->buf4_close(&Eaijk);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }

  else if(params.ref == 2) { /*** UHF ***/
    /** HC1_Wmnij():  Wmnij = + P(ij) Cje <mn||ie> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);

    global_dpd_->buf4_init(&WMNIJ, PSIF_CC3_HC1, C_irr, 2, 0, 2, 2, 0, "HC1 WMNIJ (M>N,I>J)");
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    global_dpd_->contract424(&Eijka, &CME, &WMNIJ, 3, 1, 0, 1, 0.0);
    global_dpd_->contract244(&CME, &Eaijk, &WMNIJ, 1, 0, 1, 1, 1.0);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk);
    global_dpd_->buf4_close(&WMNIJ);

    global_dpd_->buf4_init(&Wmnij, PSIF_CC3_HC1, C_irr, 12, 10, 12, 12, 0, "HC1 Wmnij (m>n,i>j)");
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    global_dpd_->contract424(&Eijka, &Cme, &Wmnij, 3, 1, 0, 1, 0.0);
    global_dpd_->contract244(&Cme, &Eaijk, &Wmnij, 1, 0, 1, 1, 1.0);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk);
    global_dpd_->buf4_close(&Wmnij);

    global_dpd_->buf4_init(&WMnIj, PSIF_CC3_HC1, C_irr, 22, 22, 22, 22, 0, "HC1 WMnIj (Mn,Ij)");
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
    global_dpd_->contract424(&Eijka, &Cme, &WMnIj, 3, 1, 0, 1, 0.0);
    global_dpd_->contract244(&CME, &Eaijk, &WMnIj, 1, 0, 1, 1, 1.0);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk);
    global_dpd_->buf4_close(&WMnIj);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }
}


void HC1_Wmbij(int i, int C_irr)
{
  dpdbuf4 W, Wmnij, I, Z, Z1, Z2, C, D;
  dpdfile2 Cme, CME;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if(params.ref == 0) { /** RHF **/
    /* HC1_Wmbij = +P(ij) Cie <mb||ej> - Cnb <mn||ij> */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");

    /** - C_n^b <Mn|Ij> -> W(Mb,Ij) **/
    global_dpd_->buf4_init(&I, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->contract424(&I, &CME, &W, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&I);

    /* <Mb||Ie> Cje -> W(Mb,Ij) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &CME, &W, 3, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&C);

    /* CIE <Mb|Ej> -> W(Mb,Ij) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&CME, &D, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, rspq, 0, 10, "HC1 WMbIj (Ij,Mb)");

    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&CME);
  }

  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);

    /** - C_N^B <MN||IJ> --> W(MB,IJ) **/
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 2, 10, 2, 0, "HC1 WMBIJ (MB,I>J)");
    global_dpd_->contract424(&Wmnij, &CME, &W, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);

    /** - C_n^b <mn||ij> **/
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 2, 10, 2, 0, "HC1 Wmbij (mb,i>j)");
    global_dpd_->contract424(&Wmnij, &Cme, &W, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);

    /** - C_n^b <Mn|Ij> **/
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");
    global_dpd_->contract424(&Wmnij, &Cme, &W, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);

    /** - C_N^B <mN|iJ> **/
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WmBiJ (mB,iJ)");
    global_dpd_->contract424(&Wmnij, &CME, &W, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);

    /* term 2: - P(ij) <mb||ei> Cje -> Wmbij */

    /** + P(IJ) C_J^E <MB||IE> -> WMBIJ **/
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 10, 0, 10, 0, 0, "Z1(MB,IJ)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &CME, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, pqsr, 10, 0, "Z2(MB,JI)");

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 10, 0, 10, 0, 0, "Z2(MB,JI)");
    global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0, 10, 2, 0, "HC1 WMBIJ (MB,I>J)");
    global_dpd_->buf4_axpy(&Z2, &W, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);

    /** - P(ij) C_j^e ( <mb||ie> ) -> WMBIJ **/

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 10, 0, 10, 0, 0, "Z1(mb,ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &Cme, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, pqsr, 10, 0, "Z2(mb,ji)");

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 10, 0, 10, 0, 0, "Z2(mb,ji)");
    global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0, 10, 2, 0, "HC1 Wmbij (mb,i>j)");
    global_dpd_->buf4_axpy(&Z2, &W, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);

    /* <Mb||Ie> Cje -> W(Mb,Ij) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &Cme, &W, 3, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&C);

    /* CIE <Mb|Ej> -> W(Mb,Ij) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&CME, &D, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /** C_J^E <mB||iE> = Cje <mB|iE> **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WmBiJ (mB,iJ)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &CME, &W, 3, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&C);

    /** -C_i^e <mB||Je> = +Cie <mB|eJ> = +Cie <mJ|eB> **/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&Cme, &D, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /* do purge before sort */
    purge_Wmbij(C_irr);

    /* do final sort to get (Ij,Mb) */

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 2, 10, 2, 0, "HC1 WMBIJ (MB,I>J)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, rspq, 2, 10, "HC1 WMBIJ (I>J,MB)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 2, 10, 2, 0, "HC1 Wmbij (mb,i>j)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, rspq, 2, 10, "HC1 Wmbij (i>j,mb)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, rspq, 0, 10, "HC1 WMbIj (Ij,Mb)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WmBiJ (mB,iJ)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, rspq, 0, 10, "HC1 WmBiJ (iJ,mB)");
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /** - C_N^B W_MNIJ --> W(MB,IJ) **/
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_AINTS, 0, 0, 2, 0, 0, 1, "A <IJ|KL>");
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 20, 2, 20, 2, 0, "HC1 WMBIJ (MB,I>J)");
    global_dpd_->contract424(&Wmnij, &CME, &W, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);

    /** - C_n^b W_mnij **/
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_AINTS, 0, 10, 12, 10, 10, 1, "A <ij|kl>");
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 30, 12, 30, 12, 0, "HC1 Wmbij (mb,i>j)");
    global_dpd_->contract424(&Wmnij, &Cme, &W, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);

    /** - C_n^b W_MnIj **/
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 24, 22, 24, 22, 0, "HC1 WMbIj (Mb,Ij)");
    global_dpd_->contract424(&Wmnij, &Cme, &W, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);

    /** - C_N^B W_mNiJ **/
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    global_dpd_->buf4_sort(&Wmnij, PSIF_CC_TMP0, qpsr, 23, 23, "A <iJ|kL>");
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_TMP0, 0, 23, 23, 23, 23, 0, "A <iJ|kL>");
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 27, 23, 27, 23, 0, "HC1 WmBiJ (mB,iJ)");
    global_dpd_->contract424(&Wmnij, &CME, &W, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Wmnij);

    /* term 2: - P(ij) <mb||ei> Cje -> Wmbij */

    /** + P(IJ) C_J^E <MB||IE> -> WMBIJ **/
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 20, 0, 20, 0, 0, "Z1(MB,IJ)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    global_dpd_->contract424(&C, &CME, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, pqsr, 20, 0, "Z2(MB,JI)");

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 20, 0, 20, 0, 0, "Z2(MB,JI)");
    global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 20, 0, 20, 2, 0, "HC1 WMBIJ (MB,I>J)");
    global_dpd_->buf4_axpy(&Z2, &W, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);

    /** - P(ij) C_j^e ( <mb||ie> ) -> WMBIJ **/

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 30, 10, 30, 10, 0, "Z1(mb,ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &Cme, &Z2, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, pqsr, 30, 10, "Z2(mb,ji)");

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 30, 10, 30, 10, 0, "Z2(mb,ji)");
    global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 30, 10, 30, 12, 0, "HC1 Wmbij (mb,i>j)");
    global_dpd_->buf4_axpy(&Z2, &W, 1.0);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);

    /* <Mb||Ie> Cje -> W(Mb,Ij) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 24, 22, 24, 22, 0, "HC1 WMbIj (Mb,Ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    global_dpd_->contract424(&C, &Cme, &W, 3, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&C);

    /* CIE <Mb|Ej> -> W(Mb,Ij) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->contract244(&CME, &D, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /** C_J^E <mB||iE> = C^J_E <mB|iE> **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 27, 23, 27, 23, 0, "HC1 WmBiJ (mB,iJ)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    global_dpd_->contract424(&C, &CME, &W, 3, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&C);

    /** -C_i^e <mB||Je> = +Cie <mB|eJ> = +Cie <mJ|eB> **/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->contract244(&Cme, &D, &W, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /** do final sort to (Ij,Mb) **/
    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 20, 2, 20, 2, 0, "HC1 WMBIJ (MB,I>J)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, rspq, 2, 20, "HC1 WMBIJ (I>J,MB)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 30, 12, 30, 12, 0, "HC1 Wmbij (mb,i>j)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, rspq, 12, 30, "HC1 Wmbij (i>j,mb)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 24, 22, 24, 22, 0, "HC1 WMbIj (Mb,Ij)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, rspq, 22, 24, "HC1 WMbIj (Ij,Mb)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 27, 23, 27, 23, 0, "HC1 WmBiJ (mB,iJ)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, rspq, 23, 27, "HC1 WmBiJ (iJ,mB)");
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }

}


void HC1_Wmbej(int i, int C_irr) {
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ;
  dpdbuf4 D, C, F, E, X, Y, W, Z;
  dpdfile2 Cme, CME;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if(params.ref == 0) { /** RHF **/
    /* Cjf <mb||ef> -> Wmbej */

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);

    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->contract424(&F, &CME, &WMbEj, 3, 1, 0, 1, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&WMbEj);


    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->contract424(&F, &CME, &Z, 3, 1, 0, -1, 0);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 10, 10, "WMbeJ"); /* (Mb,Je) */
    global_dpd_->buf4_close(&Z);

    /* - Cnb <mn||ej> -> Wmbej */

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&E, &CME, &WMbEj, 3, 0, 1, -1, 1.0);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&E, &CME, &WMbeJ, 1, 0, 1, 1, 1.0);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_close(&E);

    global_dpd_->file2_close(&CME);

    /* Sort to (ME,JB) */

    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC3_HC1, prsq, 10, 10, "HC1 WMbEj (ME,jb)");
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC3_HC1, psrq, 10, 10, "HC1 WMbeJ (Me,Jb)");
    global_dpd_->buf4_close(&WMbeJ);

  }
  else if(params.ref == 1) { /** ROHF **/


    /* + Cjf <mb||ef> -> Wmbej*/
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->contract424(&F, &CME, &WMBEJ, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->contract424(&F, &Cme, &Wmbej, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&F, &Cme, &WMbEj, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->contract424(&F, &CME, &WmBeJ, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&WmBeJ);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract244(&CME, &F, &WMbeJ, 1, 2, 1, -1, 0);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract244(&Cme, &F, &WmBEj, 1, 2, 1, -1, 0);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&F);

    /* - Cnb <mn||ej> -> Wmbej */

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->contract424(&E, &CME, &WMBEJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->contract424(&E, &Cme, &Wmbej, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&E, &Cme, &WMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->contract424(&E, &CME, &WmBeJ, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WmBeJ);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&E, &Cme, &WMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract424(&E, &CME, &WmBEj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&E);

    /* Convert to (ME,JB) for remaining terms */

    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->buf4_sort(&WMBEJ, PSIF_CC3_HC1, prsq, 10, 10, "HC1 WMBEJ (ME,JB)");
    global_dpd_->buf4_close(&WMBEJ);

    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->buf4_sort(&Wmbej, PSIF_CC3_HC1, prsq, 10, 10, "HC1 Wmbej (me,jb)");
    global_dpd_->buf4_close(&Wmbej);

    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC3_HC1, prsq, 10, 10, "HC1 WMbEj (ME,jb)");
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->buf4_sort(&WmBeJ, PSIF_CC3_HC1, prsq, 10, 10, "HC1 WmBeJ (me,JB)");
    global_dpd_->buf4_close(&WmBeJ);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC3_HC1, psrq, 10, 10, "HC1 WMbeJ (Me,Jb)");
    global_dpd_->buf4_close(&WMbeJ);

    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->buf4_sort(&WmBEj, PSIF_CC3_HC1, psrq, 10, 10, "HC1 WmBEj (mE,jB)");
    global_dpd_->buf4_close(&WmBEj);

  } /** ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* F -> Wmbej */ /* + Cjf <mb||ef> -> Wmbej*/

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract424(&F, &CME, &W, 3, 1, 0, 1, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &Cme, &W, 3, 1, 0, 1, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract424(&F, &Cme, &W, 3, 1, 0, 1, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract424(&F, &CME, &W, 3, 1, 0, 1, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract244(&CME, &F, &W, 1, 2, 1, -1, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract244(&Cme, &F, &W, 1, 2, 1, -1, 0.0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    /* - Cnb <mn||ej> -> Wmbej */

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    global_dpd_->contract424(&E, &CME, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->contract424(&E, &Cme, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->contract424(&E, &Cme, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    global_dpd_->contract424(&E, &CME, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->contract424(&E, &Cme, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->contract424(&E, &CME, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    /* Convert to (ME,JB) for remaining terms */

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, prsq, 20, 20, "HC1 WMBEJ (ME,JB)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, prsq, 30, 30, "HC1 Wmbej (me,jb)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, prsq, 20, 30, "HC1 WMbEj (ME,jb)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, prsq, 30, 20, "HC1 WmBeJ (me,JB)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, psrq, 24, 24, "HC1 WMbeJ (Me,Jb)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, psrq, 27, 27, "HC1 WmBEj (mE,jB)");
    global_dpd_->buf4_close(&W);
  }

  return;
}



void HC1_Wabei(int i, int C_irr) {
  dpdbuf4 Z, Z1, Z2, Z3;
  dpdbuf4 B, C, D, E, F, W;
  dpdfile2 Cme, CME, T1;
  char CME_lbl[32], Cme_lbl[32];
  int Gef, Gei, Gab, Ge, Gf, Gi;
  int EE, e;
  int nrows, ncols, nlinks;

  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if(params.ref == 0) { /** RHF **/


    /* Z1(Ab,Ei) <-- <Ab|Ef> * C(i,f) */
/*
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "CC3 Z(Ab,Ei)");
    dpd_buf4_init(&Z2, CC_TMP0, C_irr, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract424(&B, &CME, &Z1, 3, 1, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&Z2);
    dpd_file2_close(&CME);
*/

    // Added new B(+)/B(-) code from cchbar 12/29/09, -TDC
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
    global_dpd_->file2_init(&T1, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 11, 8, 11, 8, 0, "Z1(ei,a>=b)");
    global_dpd_->buf4_scm(&Z1, 0.0); /* needed for empty occpi or virtpi irreps */
    for(Gef=0; Gef < moinfo.nirreps; Gef++) {
      Gab = Gef; /* B is totally symmetric */
      for(Ge=0; Ge < moinfo.nirreps; Ge++) {
        Gf = Ge ^ Gef;
        Gi = Gf ^ C_irr;  /* T1 is not necessarily totally symmetric */
        Gei = Ge ^ Gi;
        if((Gei^Gab) != C_irr) printf("CC3 HC1 Error: C_irr = %d; Gef = %d; Ge = %d; Gi = %d\n", C_irr, Gef, Ge, Gi);
        B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gab]);
        Z1.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gi],Z1.params->coltot[Gab])
;
        nrows = moinfo.occpi[Gi];
        ncols = Z1.params->coltot[Gab];
        nlinks = moinfo.virtpi[Gf];
        if(nrows && ncols && nlinks) {
          for(EE=0; EE < moinfo.virtpi[Ge]; EE++) {
            e = moinfo.vir_off[Ge] + EE;
            global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
            C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,
                    B.matrix[Gef][0],ncols,0.0,Z1.matrix[Gei][0],ncols);
            global_dpd_->buf4_mat_irrep_wrt_block(&Z1,Gei,Z1.row_offset[Gei][e],moinfo.occpi[Gi]);
          }
        }
        global_dpd_->free_dpd_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gab]);
        global_dpd_->free_dpd_block(Z1.matrix[Gei], moinfo.occpi[Gi], Z1.params->coltot[Gab]);
      }
    }
    global_dpd_->buf4_close(&Z1);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&B);

    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
    global_dpd_->file2_init(&T1, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 11, 9, 11, 9, 0, "Z2(ei,a>=b)");
    global_dpd_->buf4_scm(&Z2, 0.0); /* needed for empty occpi or virtpi irreps */
    for(Gef=0; Gef < moinfo.nirreps; Gef++) {
      Gab = Gef; /* W and B are totally symmetric */
      for(Ge=0; Ge < moinfo.nirreps; Ge++) {
        Gf = Ge ^ Gef;
        Gi = Gf ^ C_irr;  /* T1 is not necessarily totally symmetric */
        Gei = Ge ^ Gi;
        if((Gei^Gab) != C_irr) printf("CC3 HC1 Error: C_irr = %d; Gef = %d; Ge = %d; Gi = %d\n", C_irr, Gef, Ge, Gi);
        B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gab]);
        Z2.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gi],Z2.params->coltot[Gab]);

        nrows = moinfo.occpi[Gi];
        ncols = Z2.params->coltot[Gab];
        nlinks = moinfo.virtpi[Gf];
        if(nrows && ncols && nlinks) {
          for(EE=0; EE < moinfo.virtpi[Ge]; EE++) {
            e = moinfo.vir_off[Ge] + EE;
            global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
            C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,
                    B.matrix[Gef][0],ncols,0.0,Z2.matrix[Gei][0],ncols);
            global_dpd_->buf4_mat_irrep_wrt_block(&Z2, Gei, Z2.row_offset[Gei][e], moinfo.occpi[Gi]);
          }
        }
        global_dpd_->free_dpd_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gab]);
        global_dpd_->free_dpd_block(Z2.matrix[Gei], moinfo.occpi[Gi], Z2.params->coltot[Gab]);
      }
    }
    global_dpd_->buf4_close(&Z2);
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&B);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 11, 5, 11, 8, 0, "Z1(ei,a>=b)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 11, 5, 11, 9, 0, "Z2(ei,a>=b)");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");
    global_dpd_->buf4_scm(&W, 0);
    global_dpd_->buf4_axpy(&Z1, &W, 1);
    global_dpd_->buf4_axpy(&Z2, &W, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "CC3 Z(Ab,Ei)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");

    /* Z1(Ab,Ei) <--  - C(M,A) * <Mb|Ei> */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&CME, &D, &Z1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&D);

    /* Z2(Ei,Ab) <-- - <mA,iE> C(m,b) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_sort(&C, PSIF_CC_TMP0, qpsr, 11, 11, "CC3 Z(Ei,Am)");
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 11, 11, 11, 0, "CC3 Z(Ei,Am)");
    global_dpd_->contract424(&Z, &CME, &Z2, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_close(&Z2);

    /* W(Ab,Ei) = Z1(Ab,Ei) + Z2(Ei,Ab) */
    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TMP0, rspq, 11, 5, "CC3 Z(Ei,Ab)", 1);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, C_irr, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");
    global_dpd_->buf4_sort(&Z2, PSIF_CC3_HC1, qpsr, 10, 5, "CC3 WAbEi (Ie,Ab)");

    global_dpd_->file2_close(&CME);
  }

  else if (params.ref == 1) { /* ROHF */
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 7, 11, 7, 11, 0, "WABEI (A>B,EI)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <ab|cd>");
    global_dpd_->contract424(&B, &CME, &Z1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 7, 11, 7, 11, 0, "Wabei (a>b,ei)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <ab|cd>");
    global_dpd_->contract424(&B, &Cme, &Z1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->contract424(&B, &Cme, &Z1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WaBeI (aB,eI)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->contract424(&B, &CME, &Z1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Z1);

    /** -CMA <MI||EB> + CMB <MI||EA> **/

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 5, 11, 5, 11, 0, "Z (AB,EI)");

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->contract244(&CME, &C, &Z1, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qprs, 5, 11, "Z (BA,EI)");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 5, 11, 7, 11, 0, "WABEI (A>B,EI)");
    global_dpd_->buf4_axpy(&Z1, &W, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 5, 11, 5, 11, 0, "Z (BA,EI)");
    global_dpd_->buf4_axpy(&Z1, &W, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&W);

    /** -Cma <mi||eb> + Cmb <mi||ea> **/

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 5, 11, 5, 11, 0, "Z (ab,ei)");

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->contract244(&Cme, &C, &Z1, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qprs, 5, 11, "Z (ba,ei)");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 5, 11, 7, 11, 0, "Wabei (a>b,ei)");
    global_dpd_->buf4_axpy(&Z1, &W, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 5, 11, 5, 11, 0, "Z (ba,ei)");
    global_dpd_->buf4_axpy(&Z1, &W, -1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_close(&W);

    /** -CMA <Mi||Eb> - Cmb <mA||iE> **/

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&CME, &D, &Z1, 0, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 5, 11, 5, 11, 0, "WAbEi (bA,Ei)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia|jb> (ia,bj)");
    global_dpd_->contract244(&Cme, &C, &Z1, 0, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TMP0, qprs, 5, 11, "WAbEi (Ab,Ei)", 1.0);
    global_dpd_->buf4_close(&Z1);

    /** **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WaBeI (aB,eI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&Cme, &D, &Z1, 0, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 5, 11, 5, 11, 0, "WaBeI (Ba,eI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia|jb> (ia,bj)");
    global_dpd_->contract244(&CME, &C, &Z1, 0, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TMP0, qprs, 5, 11, "WaBeI (aB,eI)", 1.0);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);

    /* final sort to (EI,AB) */

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 5, 11, 7, 11, 0, "WABEI (A>B,EI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP0, rspq, 11, 7, "WABEI (EI,A>B)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 5, 11, 7, 11, 0, "Wabei (a>b,ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP0, rspq, 11, 7, "Wabei (ei,a>b)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WAbEi (Ab,Ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP0, rspq, 11, 5, "WAbEi (Ei,Ab)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WaBeI (aB,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP0, rspq, 11, 5, "WaBeI (eI,aB)");
    global_dpd_->buf4_close(&W);

    purge_Wabei(C_irr);

    /* final sort to Wabei (IE,AB) */
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 11, 7, 11, 7, 0, "WABEI (EI,A>B)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, qprs, 10, 7, "HC1 WABEI (IE,A>B)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 11, 7, 11, 7, 0, "Wabei (ei,a>b)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, qprs, 10, 7, "HC1 Wabei (ie,a>b)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, qprs, 10, 5, "HC1 WAbEi (iE,Ab)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 11, 5, 11, 5, 0, "WaBeI (eI,aB)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, qprs, 10, 5, "HC1 WaBeI (Ie,aB)");
    global_dpd_->buf4_close(&W);
  }
  else if (params.ref == 2) { /* UHF */
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /* term 1, Cif <ab||ef> */

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 7, 21, 7, 21, 0, "WABEI (A>B,EI)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
    global_dpd_->contract424(&B, &CME, &Z1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 17, 31, 17, 31, 0, "Wabei (a>b,ei)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 17, 15, 15, 15, 1, "B <ab|cd>");
    global_dpd_->contract424(&B, &Cme, &Z1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 28, 26, 28, 26, 0, "WAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    global_dpd_->contract424(&B, &Cme, &Z1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&B);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 24, 28, 24, 28, 0, "Z(Ie,Ba)");
    global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    global_dpd_->contract244(&CME, &B, &Z, 1, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&B);
    /** Z(Ie,Ba) --> W'(aB,eI) **/
    /* srqp seems to have a bug
     dpd_buf4_sort(&Z, CC_TMP0, srqp, 29, 25, "WaBeI (aB,eI)");
    */
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, rspq, 28, 24, "WaBeI (Ba,Ie) 1");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, C_irr, 28, 24, 28, 24, 0, "WaBeI (Ba,Ie)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 29, 25, "WaBeI (aB,eI)");
    global_dpd_->buf4_close(&Z);

    /** UHF term 2 **/
    /* -CMA <MI||EB> + CMB <MI||EA> **/

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 5, 21, 5, 21, 0, "Z (AB,EI)");

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    global_dpd_->contract244(&CME, &C, &Z1, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qprs, 5, 21, "Z (BA,EI)");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 5, 21, 7, 21, 0, "WABEI (A>B,EI)");
    global_dpd_->buf4_axpy(&Z1, &W, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 5, 21, 5, 21, 0, "Z (BA,EI)");
    global_dpd_->buf4_axpy(&Z1, &W, -1.0);
    global_dpd_->buf4_close(&Z1);

    /** -Cma <mi||eb> + Cmb <mi||ea> **/

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 15, 31, 15, 31, 0, "Z (ab,ei)");

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->contract244(&Cme, &C, &Z1, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qprs, 15, 31, "Z (ba,ei)");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 15, 31, 17, 31, 0, "Wabei (a>b,ei)");
    global_dpd_->buf4_axpy(&Z1, &W, 1.0);
    global_dpd_->buf4_close(&Z1);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 15, 31, 15, 31, 0, "Z (ba,ei)");
    global_dpd_->buf4_axpy(&Z1, &W, -1.0);
    global_dpd_->buf4_close(&Z1);

    /** -CMA <Mi||Eb> - Cmb <mA|iE> **/

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 28, 26, 28, 26, 0, "WAbEi (Ab,Ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->contract244(&CME, &D, &Z1, 0, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 29, 26, 29, 26, 0, "WAbEi (bA,Ei)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 26, 27, 26, 0, "C <Ai|Bj> (iA,Bj)");
    global_dpd_->contract244(&Cme, &C, &Z1, 0, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TMP0, qprs, 28, 26, "WAbEi (Ab,Ei)", 1.0);
    global_dpd_->buf4_close(&Z1);

    /**  **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, C_irr, 29, 25, 29, 25, 0, "WaBeI (aB,eI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->contract244(&Cme, &D, &Z1, 0, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, C_irr, 28, 25, 28, 25, 0, "WaBeI (Ba,eI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 25, 24, 25, 0, "C <Ia|Jb> (Ia,bJ)");
    global_dpd_->contract244(&CME, &C, &Z1, 0, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_TMP0, qprs, 29, 25, "WaBeI (aB,eI)", 1.0);
    global_dpd_->buf4_close(&Z1);

    /* final sort and storage */
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 7, 21, 7, 21, 0, "WABEI (A>B,EI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP0, rspq, 21, 7, "WABEI (EI,A>B)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 21, 7, 21, 7, 0, "WABEI (EI,A>B)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, qprs, 20, 7, "HC1 WABEI (IE,B>A)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 20, 7, 20, 7, 0, "HC1 WABEI (IE,B>A)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);

    /* final sort and storage */
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 17, 31, 17, 31, 0, "Wabei (a>b,ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP0, rspq, 31, 17, "Wabei (ei,a>b)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 31, 17, 31, 17, 0, "Wabei (ei,a>b)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, qprs, 30, 17, "HC1 Wabei (ie,b>a)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC3_HC1, C_irr, 30, 17, 30, 17, 0, "HC1 Wabei (ie,b>a)");
    global_dpd_->buf4_scm(&W, -1.0);
    global_dpd_->buf4_close(&W);

    /* final sort and storage */
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 28, 26, 28, 26, 0, "WAbEi (Ab,Ei)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP0, rspq, 26, 28, "WAbEi (Ei,Ab)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 26, 28, 26, 28, 0, "WAbEi (Ei,Ab)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, qpsr, 27, 29, "HC1 WAbEi (iE,bA)");
    global_dpd_->buf4_close(&W);

    /* final sort and storage */
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 29, 25, 29, 25, 0, "WaBeI (aB,eI)");
    global_dpd_->buf4_sort(&W, PSIF_CC_TMP0, rspq, 25, 29, "WaBeI (eI,aB)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, C_irr, 25, 29, 25, 29, 0, "WaBeI (eI,aB)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HC1, qpsr, 24, 28, "HC1 WaBeI (Ie,Ba)");
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
  }
}

void purge_HC1(int C_irr) {
  dpdfile2 FAE, Fmi, FME, Fme;
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  /* Purge FME matrix elements */
  global_dpd_->file2_init(&FME, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 FME");
  global_dpd_->file2_mat_init(&FME);
  global_dpd_->file2_mat_rd(&FME);
  for(h=0; h < nirreps; h++) {
    for(m=0; m<occpi[h]; m++)
      for(e=(virtpi[h]-openpi[h]); e<virtpi[h]; e++)
	FME.matrix[h][m][e] = 0.0;
  }
  global_dpd_->file2_mat_wrt(&FME);
  global_dpd_->file2_mat_close(&FME);
  global_dpd_->file2_close(&FME);

  /* Purge Fme matrix elements */
  global_dpd_->file2_init(&Fme, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 Fme");
  global_dpd_->file2_mat_init(&Fme);
  global_dpd_->file2_mat_rd(&Fme);
  for(h=0; h < nirreps; h++) {
    for(e=0; e<virtpi[h]; e++)
      for(m=(occpi[h]-openpi[h]); m<occpi[h]; m++)
	Fme.matrix[h][m][e] = 0.0;
  }
  global_dpd_->file2_mat_wrt(&Fme);
  global_dpd_->file2_mat_close(&Fme);
  global_dpd_->file2_close(&Fme);

  /* Purge Wmnij matrix elements */
  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 2, 2,"HC1 Wmnij (m>n,i>j)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn < W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ij=0; ij < W.params->coltot[h]; ij++) {
	i = W.params->colorb[h][ij][0];
	j = W.params->colorb[h][ij][1];
	isym = W.params->rsym[i];
	jsym = W.params->ssym[j];
	I = i - occ_off[isym];
	J = j - occ_off[jsym];
	if ((I >= (occpi[isym] - openpi[isym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) ||
	    (M >= (occpi[msym] - openpi[msym])) ||
	    (N >= (occpi[nsym] - openpi[nsym])) )
	  W.matrix[h][mn][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 0, 0,"HC1 WMnIj (Mn,Ij)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn < W.params->rowtot[h]; mn++) {
      n = W.params->roworb[h][mn][1];
      nsym = W.params->qsym[n];
      N = n - occ_off[nsym];
      for(ij=0; ij < W.params->coltot[h]; ij++) {
	j = W.params->colorb[h][ij][1];
	jsym = W.params->ssym[j];
	J = j - occ_off[jsym];
	if ((J >= (occpi[jsym] - openpi[jsym])) ||
	    (N >= (occpi[nsym] - openpi[nsym])) )
	  W.matrix[h][mn][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);


  /* Purge Wmbej matrix elements */
  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 10, 10,"HC1 WMBEJ (ME,JB)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      e = W.params->roworb[h][me][1];
      esym = W.params->qsym[e];
      E = e - vir_off[esym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	b = W.params->colorb[h][jb][1];
	bsym = W.params->ssym[b];
	B = b - vir_off[bsym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);


  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 10, 10,"HC1 Wmbej (me,jb)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      m = W.params->roworb[h][me][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	j = W.params->colorb[h][jb][0];
	jsym = W.params->rsym[j];
	J = j - occ_off[jsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);


  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 10, 10,"HC1 WMbEj (ME,jb)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      e = W.params->roworb[h][me][1];
      esym = W.params->qsym[e];
      E = e - vir_off[esym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	j = W.params->colorb[h][jb][0];
	jsym = W.params->rsym[j];
	J = j - occ_off[jsym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);


  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 10, 10,"HC1 WmBeJ (me,JB)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      m = W.params->roworb[h][me][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	b = W.params->colorb[h][jb][1];
	bsym = W.params->ssym[b];
	B = b - vir_off[bsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);


  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 10, 10,"HC1 WmBEj (mE,jB)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(me=0; me < W.params->rowtot[h]; me++) {
      m = W.params->roworb[h][me][0];
      e = W.params->roworb[h][me][1];
      msym = W.params->psym[m];
      esym = W.params->qsym[e];
      M = m - occ_off[msym];
      E = e - vir_off[esym];
      for(jb=0; jb < W.params->coltot[h]; jb++) {
	j = W.params->colorb[h][jb][0];
	b = W.params->colorb[h][jb][1];
	jsym = W.params->rsym[j];
	bsym = W.params->ssym[b];
	J = j - occ_off[jsym];
	B = b - vir_off[bsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (E >= (virtpi[esym] - openpi[esym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][me][jb] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);


  /* Purge Wamef matrix elements */
  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 11, 7,"HC1 WAMEF (AM,E>F)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      a = W.params->roworb[h][ma][0];
      asym = W.params->psym[a];
      A = a - vir_off[asym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	e = W.params->colorb[h][ef][0];
	f = W.params->colorb[h][ef][1];
	esym = W.params->rsym[e];
	fsym = W.params->ssym[f];
	E = e - vir_off[esym];
	F = f - vir_off[fsym];
	if ((A >= (virtpi[asym] - openpi[asym])) ||
	    (E >= (virtpi[esym] - openpi[esym])) ||
	    (F >= (virtpi[fsym] - openpi[fsym])) )
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 11, 7,"HC1 Wamef (am,e>f)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      m = W.params->roworb[h][ma][1];
      msym = W.params->qsym[m];
      M = m - occ_off[msym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	if (M >=  (occpi[msym] - openpi[msym]))
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 11, 5,"HC1 WAmEf (Am,Ef)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      a = W.params->roworb[h][ma][0];
      m = W.params->roworb[h][ma][1];
      asym = W.params->psym[a];
      msym = W.params->qsym[m];
      M = m - occ_off[msym];
      A = a - vir_off[asym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	e = W.params->colorb[h][ef][0];
	esym = W.params->rsym[e];
	E = e - vir_off[esym];
	if ((A >= (virtpi[asym] - openpi[asym])) ||
	    (M >=  (occpi[msym] - openpi[msym])) ||
	    (E >= (virtpi[esym] - openpi[esym])) )
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 11, 5,"HC1 WaMeF (aM,eF)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	f = W.params->colorb[h][ef][1];
	fsym = W.params->ssym[f];
	F = f - vir_off[fsym];
	if (F >= (virtpi[fsym] - openpi[fsym]))
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);



  return;
}

/* Purge Wmnie matrix elements */
void purge_Wmnie(int C_irr) {
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 0, 11,"HC1 WMnIe (Mn,eI)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      n = W.params->roworb[h][mn][1];
      nsym = W.params->qsym[n];
      N = n - occ_off[nsym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
	if (N >= (occpi[nsym] - openpi[nsym]))
	  W.matrix[h][mn][ei] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 2, 11, "HC1 WMNIE (M>N,EI)");
  for(h=0; h < W.params->nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        e = W.params->colorb[h][ei][0];
        esym = W.params->rsym[e];
        E = e - vir_off[esym];
        if (E >= (virtpi[esym] - openpi[esym]))
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 2, 11,"HC1 Wmnie (m>n,ei)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      n = W.params->roworb[h][mn][1];
      msym = W.params->psym[m];
      nsym = W.params->qsym[n];
      M = m - occ_off[msym];
      N = n - occ_off[nsym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        i = W.params->colorb[h][ei][1];
        isym = W.params->ssym[i];
        I = i - occ_off[isym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
          (N >= (occpi[nsym] - openpi[nsym])) ||
          (I >= (occpi[isym] - openpi[isym])) )
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 0, 11,"HC1 WmNiE (mN,Ei)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      m = W.params->roworb[h][mn][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        e = W.params->colorb[h][ei][0];
        i = W.params->colorb[h][ei][1];
        esym = W.params->rsym[e];
        isym = W.params->ssym[i];
        E = e - vir_off[esym];
        I = i - occ_off[isym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
            (E >= (virtpi[esym] - openpi[esym])) ||
            (I >= (occpi[isym] - openpi[isym])) )
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);
  return;
}

/* Purge WMBIJ matrix elements */
void purge_Wmbij(int C_irr) {
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 10, 2,"HC1 WMBIJ (MB,I>J)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      b = W.params->roworb[h][mb][1];
      bsym = W.params->qsym[b];
      B = b - vir_off[bsym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	if (B >= (virtpi[bsym] - openpi[bsym]))
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 10, 2,"HC1 Wmbij (mb,i>j)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      m = W.params->roworb[h][mb][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	i = W.params->colorb[h][ij][0];
	j = W.params->colorb[h][ij][1];
	isym = W.params->rsym[i];
	jsym = W.params->ssym[j];
	I = i - occ_off[isym];
	J = j - occ_off[jsym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (I >= (occpi[isym] - openpi[isym])) ||
	    (J >= (occpi[jsym] - openpi[jsym])) )
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0,"HC1 WMbIj (Mb,Ij)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	j = W.params->colorb[h][ij][1];
	jsym = W.params->ssym[j];
	J = j - occ_off[jsym];
	if (J >= (occpi[jsym] - openpi[jsym]))
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC3_HC1, C_irr, 10, 0,"HC1 WmBiJ (mB,iJ)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      m = W.params->roworb[h][mb][0];
      b = W.params->roworb[h][mb][1];
      msym = W.params->psym[m];
      bsym = W.params->qsym[b];
      M = m - occ_off[msym];
      B = b - vir_off[bsym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	i = W.params->colorb[h][ij][0];
	isym = W.params->rsym[i];
	I = i - occ_off[isym];
	if ((M >= (occpi[msym] - openpi[msym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) ||
	    (I >= (occpi[isym] - openpi[isym])) )
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);
}


/* Purge Wabei matrix elements */
void purge_Wabei(int C_irr) {

  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  global_dpd_->file4_init(&W, PSIF_CC_TMP0, C_irr, 11, 7,"WABEI (EI,A>B)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      e = W.params->roworb[h][ei][0];
      esym = W.params->psym[e];
      E = e - vir_off[esym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	a = W.params->colorb[h][ab][0];
	b = W.params->colorb[h][ab][1];
	asym = W.params->rsym[a];
	bsym = W.params->ssym[b];
	A = a - vir_off[asym];
	B = b - vir_off[bsym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (A >= (virtpi[asym] - openpi[asym])) ||
	    (B >= (virtpi[bsym] - openpi[bsym])) )
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC_TMP0, C_irr, 11, 7,"Wabei (ei,a>b)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      i = W.params->roworb[h][ei][1];
      isym = W.params->qsym[i];
      I = i - occ_off[isym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	if (I >= (occpi[isym] - openpi[isym]))
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC_TMP0, C_irr, 11, 5,"WAbEi (Ei,Ab)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      e = W.params->roworb[h][ei][0];
      i = W.params->roworb[h][ei][1];
      esym = W.params->psym[e];
      isym = W.params->qsym[i];
      E = e - vir_off[esym];
      I = i - occ_off[isym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	a = W.params->colorb[h][ab][0];
	asym = W.params->rsym[a];
	bsym = W.params->ssym[b];
	A = a - vir_off[asym];
	if ((E >= (virtpi[esym] - openpi[esym])) ||
	    (I >= (occpi[isym] - openpi[isym])) ||
	    (A >= (virtpi[asym] - openpi[asym])) )
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);

  global_dpd_->file4_init(&W, PSIF_CC_TMP0, C_irr, 11, 5,"WaBeI (eI,aB)");
  for(h=0; h < nirreps; h++) {
    global_dpd_->file4_mat_irrep_init(&W, h);
    global_dpd_->file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	b = W.params->colorb[h][ab][1];
	bsym = W.params->ssym[b];
	B = b - vir_off[bsym];
	if (B >= (virtpi[bsym] - openpi[bsym]))
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    global_dpd_->file4_mat_irrep_wrt(&W, h);
    global_dpd_->file4_mat_irrep_close(&W, h);
  }
  global_dpd_->file4_close(&W);
}

}} // namespace psi::cceom
