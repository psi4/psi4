/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
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

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);

    dpd_file2_init(&FME, CC3_HC1, C_irr, 0, 1, "HC1 FME");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_dot13(&CME, &D, &FME, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_file2_close(&FME);

    dpd_file2_close(&CME);
  }

  else if(params.ref == 1) { /** ROHF **/
    /* HC1_F()  Fme = +C_n^f <mn||ef> */

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    dpd_file2_init(&FME, CC3_HC1, C_irr, 0, 1, "HC1 FME");
    dpd_file2_init(&Fme, CC3_HC1, C_irr, 0, 1, "HC1 Fme");
  
    dpd_buf4_init(&D_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_dot13(&CME, &D_anti, &FME, 0, 0, 1.0, 0.0);
    dpd_dot13(&Cme, &D,      &FME, 0, 0, 1.0, 1.0);
    dpd_dot13(&Cme, &D_anti, &Fme, 0, 0, 1.0, 0.0);
    dpd_dot13(&CME, &D,      &Fme, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&D_anti);
    dpd_buf4_close(&D);

    dpd_file2_close(&FME);
    dpd_file2_close(&Fme);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/
    /* HC1_F()  Fme = +C_n^f <mn||ef> */

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    dpd_file2_init(&FME, CC3_HC1, C_irr, 0, 1, "HC1 FME");
    dpd_file2_init(&Fme, CC3_HC1, C_irr, 2, 3, "HC1 Fme");

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_contract422(&D, &CME, &FME, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_contract422(&D, &Cme, &FME, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_contract422(&D, &Cme, &Fme, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_contract422(&D, &CME, &Fme, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&D);

    dpd_file2_close(&FME);
    dpd_file2_close(&Fme);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
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

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);

    dpd_buf4_init(&W, CC3_HC1, C_irr, 11, 5, 11, 5, 0, "HC1 WAmEf (Am,Ef)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&CME, &D, &W, 0, 0, 0, -1, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HC1, qprs, 10, 5, "HC1 WAmEf (mA,Ef)");
    dpd_buf4_close(&W);

    dpd_file2_close(&CME);
  }

  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    /* C(N,A) <NM||EF> --> W(AM,E>F) */
    dpd_buf4_init(&WAMEF, CC3_HC1, C_irr, 11, 7, 11, 7, 0, "HC1 WAMEF (AM,E>F)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&CME, &D, &WAMEF, 0, 0, 0, -1, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WAMEF);  

    /* C(n,a) <nm||ef> --> W(am,e>f) */
    dpd_buf4_init(&Wamef, CC3_HC1, C_irr, 11, 7, 11, 7, 0, "HC1 Wamef (am,e>f)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&Cme, &D, &Wamef, 0, 0, 0, -1, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Wamef); 

    /* C(N,A) <Nm|Ef> --> W(Am,Ef) */
    dpd_buf4_init(&WAmEf, CC3_HC1, C_irr, 11, 5, 11, 5, 0, "HC1 WAmEf (Am,Ef)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&CME, &D, &WAmEf, 0, 0, 0, -1, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WAmEf);

    /* C(n,a) <nM|eF> --> W(aM,eF) */
    dpd_buf4_init(&WaMeF, CC3_HC1, C_irr, 11, 5, 11, 5, 0, "HC1 WaMeF (aM,eF)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&Cme, &D, &WaMeF, 0, 0, 0, -1, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WaMeF);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);

  } /** ROHF **/
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /* T(N,A) <NM||EF> --> W(AM,E>F) */
    dpd_buf4_init(&WAMEF, CC3_HC1, C_irr, 21, 7, 21, 7, 0, "HC1 WAMEF (AM,E>F)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_contract244(&CME, &D, &WAMEF, 0, 0, 0, -1, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WAMEF);  

    /* T(n,a) <nm||ef> --> W(am,e>f) */
    dpd_buf4_init(&Wamef, CC3_HC1, C_irr, 31, 17, 31, 17, 0, "HC1 Wamef (am,e>f)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&Cme, &D, &Wamef, 0, 0, 0, -1, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Wamef); 

    /* T(N,A) <Nm|Ef> --> W(Am,Ef) */
    dpd_buf4_init(&WAmEf, CC3_HC1, C_irr, 26, 28, 26, 28, 0, "HC1 WAmEf (Am,Ef)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract244(&CME, &D, &WAmEf, 0, 0, 0, -1, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WAmEf);

    /* T(n,a) <nM|eF> --> W(aM,eF) */
    dpd_buf4_init(&WaMeF, CC3_HC1, C_irr, 25, 29, 25, 29, 0, "HC1 WaMeF (aM,eF)");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract244(&Cme, &D, &WaMeF, 0, 0, 0, -1, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WaMeF);

    /* Do final sort */
    dpd_buf4_init(&WAMEF, CC3_HC1, C_irr, 21, 7, 21, 7, 0, "HC1 WAMEF (AM,E>F)");
    dpd_buf4_sort(&WAMEF, CC3_HC1, qprs, 20, 7, "HC1 WAMEF (MA,F>E)");
    dpd_buf4_close(&WAMEF);  
    dpd_buf4_init(&WAMEF, CC3_HC1, C_irr, 20, 7, 20, 7, 0, "HC1 WAMEF (MA,F>E)");
    dpd_buf4_scm(&WAMEF, -1.0);
    dpd_buf4_close(&WAMEF);  

    dpd_buf4_init(&Wamef, CC3_HC1, C_irr, 31, 17, 31, 17, 0, "HC1 Wamef (am,e>f)");
    dpd_buf4_sort(&Wamef, CC3_HC1, qprs, 30, 17, "HC1 Wamef (ma,f>e)");
    dpd_buf4_close(&Wamef); 
    dpd_buf4_init(&Wamef, CC3_HC1, C_irr, 30, 17, 30, 17, 0, "HC1 Wamef (ma,f>e)");
    dpd_buf4_scm(&Wamef, -1.0);
    dpd_buf4_close(&Wamef);  

    dpd_buf4_init(&WAmEf, CC3_HC1, C_irr, 26, 28, 26, 28, 0, "HC1 WAmEf (Am,Ef)");
    dpd_buf4_sort(&WAmEf, CC3_HC1, qpsr, 27, 29, "HC1 WAmEf (mA,fE)");
    dpd_buf4_close(&WAmEf);

    dpd_buf4_init(&WaMeF, CC3_HC1, C_irr, 25, 29, 25, 29, 0, "HC1 WaMeF (aM,eF)");
    dpd_buf4_sort(&WaMeF, CC3_HC1, qpsr, 24, 28, "HC1 WaMeF (Ma,Fe)");
    dpd_buf4_close(&WaMeF);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
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

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);

    /* C(I,F) * D(Mn,Fe) --> W(Mn,Ie) */
    dpd_buf4_init(&WMnIe, CC3_HC1, C_irr, 0, 10, 0, 10, 0, "HC1 WMnIe (Mn,Ie)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&CME, &D, &WMnIe, 1, 2, 1, 1, 0.0);
    dpd_file2_close(&CME);
    dpd_buf4_close(&D);
    /* W(Mn,Ie) --> W(Mn,eI) */
    /* dpd_buf4_sort(&WMnIe, CC3_HC1, pqsr, 0, 11, "HC1 WMnIe (Mn,eI)"); */
    dpd_buf4_close(&WMnIe);
  }

  else if(params.ref == 1) { /** ROHF **/
    /* HC1_Wmnie():  Wmnie = + Cif <mn||fe> */

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    /* D(M>N,EF) * C(I,F) --> W(M>N,EI) */
    dpd_buf4_init(&WMNIE, CC3_HC1, C_irr, 2, 11, 2, 11, 0, "HC1 WMNIE (M>N,EI)");
    dpd_buf4_init(&D_a, CC_DINTS, 0, 2, 5, 2, 5,0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D_a,&CME,&WMNIE, 3, 1, 0, -1, 0);
    dpd_buf4_close(&D_a);
    dpd_buf4_close(&WMNIE);

    /* D(m>n,ef) * C(i,f) --> W(m>n,ei) */
    dpd_buf4_init(&Wmnie, CC3_HC1, C_irr, 2, 11, 2, 11, 0, "HC1 Wmnie (m>n,ei)");
    dpd_buf4_init(&D_a, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D_a, &Cme, &Wmnie, 3, 1, 0, -1, 0);
    dpd_buf4_close(&D_a);
    dpd_buf4_close(&Wmnie);

    /* D(Mn,Fe) * C(I,F) --> W(Mn,Ie) */
    dpd_buf4_init(&WMnIe, CC_TMP0, C_irr, 0, 10, 0, 10, 0, "HC1 WMnIe (Mn,Ie)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&CME, &D, &WMnIe, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    /* W(Mn,Ie) --> W(Mn,eI) */
    dpd_buf4_sort(&WMnIe, CC3_HC1, pqsr, 0, 11, "HC1 WMnIe (Mn,eI)");
    dpd_buf4_close(&WMnIe);

    /* D(mN,fE) * C(i,f) --> W(mN,iE) */
    dpd_buf4_init(&WmNiE, CC_TMP1, C_irr, 0, 10, 0, 10, 0, "HC1 WmNiE (mN,iE)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&Cme,&D,&WmNiE, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    /* W(mN,iE) --> W(mN,Ei) */
    dpd_buf4_sort(&WmNiE, CC3_HC1, pqsr, 0, 11, "HC1 WmNiE (mN,Ei)");
    dpd_buf4_close(&WmNiE);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);

    purge_Wmnie(C_irr); /* before sorting here */

    /* also put "normal" sorted versions in CC_HBAR */
    dpd_buf4_init(&WMNIE, CC3_HC1, C_irr, 2, 11, 2, 11, 0, "HC1 WMNIE (M>N,EI)");
    dpd_buf4_sort(&WMNIE, CC3_HC1, pqsr, 2, 10, "HC1 WMNIE (M>N,IE)");
    dpd_buf4_close(&WMNIE);
    dpd_buf4_init(&Wmnie, CC3_HC1, C_irr, 2, 11, 2, 11, 0, "HC1 Wmnie (m>n,ei)");
    dpd_buf4_sort(&Wmnie, CC3_HC1, pqsr, 2, 10, "HC1 Wmnie (m>n,ie)");
    dpd_buf4_close(&Wmnie);
    dpd_buf4_init(&WMnIe, CC3_HC1, C_irr, 0, 11, 0, 11, 0, "HC1 WMnIe (Mn,eI)");
    dpd_buf4_sort(&WMnIe, CC3_HC1, pqsr, 0, 10, "HC1 WMnIe (Mn,Ie)");
    dpd_buf4_close(&WMnIe);
    dpd_buf4_init(&WmNiE, CC3_HC1, C_irr, 0, 11, 0, 11, 0, "HC1 WmNiE (mN,Ei)");
    dpd_buf4_sort(&WmNiE, CC3_HC1, pqsr, 0, 10, "HC1 WmNiE (mN,iE)");
    dpd_buf4_close(&WmNiE);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /* <M>N||EF> T(I,F) --> W(M>N,EI) */
    dpd_buf4_init(&W, CC3_HC1, C_irr, 2, 21, 2, 21, 0, "HC1 WMNIE (M>N,EI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_contract424(&D, &CME, &W, 3, 1, 0, -1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* <m>n||ef> T(i,f) --> W(m>n,ei) */
    dpd_buf4_init(&W, CC3_HC1, C_irr, 12, 31, 12, 31, 0, "HC1 Wmnie (m>n,ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D, &Cme, &W, 3, 1, 0, -1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* Z(nM,eI) = <nM|eF> T(I,F) */
    dpd_buf4_init(&Z, CC_TMP1, C_irr, 23, 25, 23, 25, 0, "Z(nM,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract424(&D, &CME, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /* Z(nM,eI) --> W(Mn,eI) */
    dpd_buf4_sort(&Z, CC3_HC1, qprs, 22, 25, "HC1 WMnIe (Mn,eI)");
    dpd_buf4_close(&Z);

    /* Z(Nm,Ei) = <Nm|Ef> T(i,f) */
    dpd_buf4_init(&Z, CC_TMP1, C_irr, 22, 26, 22, 26, 0, "Z(Nm,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &Cme, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /* Z(Nm,Ei) --> W(mN,Ei) */
    dpd_buf4_sort(&Z, CC3_HC1, qprs, 23, 26, "HC1 WmNiE (mN,Ei)");
    dpd_buf4_close(&Z);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);

    /* also put "normal" sorted versions in CC3_HC1 */
    dpd_buf4_init(&WMNIE, CC3_HC1, C_irr, 2, 21, 2, 21, 0, "HC1 WMNIE (M>N,EI)");
    dpd_buf4_sort(&WMNIE, CC3_HC1, pqsr, 2, 20, "HC1 WMNIE (M>N,IE)");
    dpd_buf4_close(&WMNIE);
    dpd_buf4_init(&Wmnie, CC3_HC1, C_irr, 12, 31, 12, 31, 0, "HC1 Wmnie (m>n,ei)");
    dpd_buf4_sort(&Wmnie, CC3_HC1, pqsr, 12, 30, "HC1 Wmnie (m>n,ie)");
    dpd_buf4_close(&Wmnie);
    dpd_buf4_init(&WMnIe, CC3_HC1, C_irr, 22, 25, 22, 25, 0, "HC1 WMnIe (Mn,eI)");
    dpd_buf4_sort(&WMnIe, CC3_HC1, pqsr, 22, 24, "HC1 WMnIe (Mn,Ie)");
    dpd_buf4_close(&WMnIe);
    dpd_buf4_init(&WmNiE, CC3_HC1, C_irr, 23, 26, 23, 26, 0, "HC1 WmNiE (mN,Ei)");
    dpd_buf4_sort(&WmNiE, CC3_HC1, pqsr, 23, 27, "HC1 WmNiE (mN,iE)");
    dpd_buf4_close(&WmNiE);
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

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);

    dpd_buf4_init(&WMnIj, CC3_HC1, C_irr, 0, 0, 0, 0, 0, "HC1 WMnIj (Mn,Ij)");
    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_contract244(&CME, &Eaijk, &WMnIj, 1, 0, 1, 1, 0.0);
    dpd_buf4_close(&Eaijk);

    dpd_buf4_init(&Eijka, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_contract424(&Eijka, &CME, &WMnIj, 3, 1, 0, 1, 1.0);
    dpd_buf4_close(&Eijka);
    dpd_buf4_close(&WMnIj);

    dpd_file2_close(&CME);
  }

  else if(params.ref == 1) { /** ROHF **/  
    /** HC1_Wmnij():  Wmnij = + P(ij) Cje <mn||ie> */

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    dpd_buf4_init(&Eijka_anti, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_init(&Eijka, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&Eaijk_anti, CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

    dpd_buf4_init(&WMNIJ, CC3_HC1, C_irr, 2, 0, 2, 2, 0, "HC1 WMNIJ (M>N,I>J)");
    dpd_contract424(&Eijka_anti, &CME, &WMNIJ, 3, 1, 0, 1, 0.0);
    dpd_contract244(&CME, &Eaijk_anti, &WMNIJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMNIJ);

    dpd_buf4_init(&Wmnij, CC3_HC1, C_irr, 2, 0, 2, 2, 0, "HC1 Wmnij (m>n,i>j)");
    dpd_contract424(&Eijka_anti, &Cme, &Wmnij, 3, 1, 0, 1, 0.0);
    dpd_contract244(&Cme, &Eaijk_anti, &Wmnij, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Wmnij);

    dpd_buf4_init(&WMnIj, CC3_HC1, C_irr, 0, 0, 0, 0, 0, "HC1 WMnIj (Mn,Ij)");
    dpd_contract424(&Eijka, &Cme, &WMnIj, 3, 1, 0, 1, 0.0);
    dpd_contract244(&CME, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMnIj);

    dpd_buf4_close(&Eijka_anti);
    dpd_buf4_close(&Eijka);
    dpd_buf4_close(&Eaijk_anti);
    dpd_buf4_close(&Eaijk);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }

  else if(params.ref == 2) { /*** UHF ***/
    /** HC1_Wmnij():  Wmnij = + P(ij) Cje <mn||ie> */

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    dpd_buf4_init(&WMNIJ, CC3_HC1, C_irr, 2, 0, 2, 2, 0, "HC1 WMNIJ (M>N,I>J)");
    dpd_buf4_init(&Eijka, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    dpd_contract424(&Eijka, &CME, &WMNIJ, 3, 1, 0, 1, 0.0);
    dpd_contract244(&CME, &Eaijk, &WMNIJ, 1, 0, 1, 1, 1.0);
    dpd_buf4_close(&Eijka);
    dpd_buf4_close(&Eaijk);
    dpd_buf4_close(&WMNIJ);

    dpd_buf4_init(&Wmnij, CC3_HC1, C_irr, 12, 10, 12, 12, 0, "HC1 Wmnij (m>n,i>j)");
    dpd_buf4_init(&Eijka, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    dpd_contract424(&Eijka, &Cme, &Wmnij, 3, 1, 0, 1, 0.0);
    dpd_contract244(&Cme, &Eaijk, &Wmnij, 1, 0, 1, 1, 1.0);
    dpd_buf4_close(&Eijka);
    dpd_buf4_close(&Eaijk);
    dpd_buf4_close(&Wmnij);

    dpd_buf4_init(&WMnIj, CC3_HC1, C_irr, 22, 22, 22, 22, 0, "HC1 WMnIj (Mn,Ij)");
    dpd_buf4_init(&Eijka, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_init(&Eaijk, CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
    dpd_contract424(&Eijka, &Cme, &WMnIj, 3, 1, 0, 1, 0.0);
    dpd_contract244(&CME, &Eaijk, &WMnIj, 1, 0, 1, 1, 1.0);
    dpd_buf4_close(&Eijka);
    dpd_buf4_close(&Eaijk);
    dpd_buf4_close(&WMnIj);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
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

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");

    /** - C_n^b <Mn|Ij> -> W(Mb,Ij) **/
    dpd_buf4_init(&I, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_contract424(&I, &CME, &W, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&I);

    /* <Mb||Ie> Cje -> W(Mb,Ij) */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_contract424(&C, &CME, &W, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&C);

    /* CIE <Mb|Ej> -> W(Mb,Ij) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&CME, &D, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&D);

    dpd_buf4_sort(&W, CC3_HC1, rspq, 0, 10, "HC1 WMbIj (Ij,Mb)");

    dpd_buf4_close(&W);
    dpd_file2_close(&CME);
  }

  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    /** - C_N^B <MN||IJ> --> W(MB,IJ) **/
    dpd_buf4_init(&Wmnij, CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 2, 10, 2, 0, "HC1 WMBIJ (MB,I>J)");
    dpd_contract424(&Wmnij, &CME, &W, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);

    /** - C_n^b <mn||ij> **/
    dpd_buf4_init(&Wmnij, CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 2, 10, 2, 0, "HC1 Wmbij (mb,i>j)");
    dpd_contract424(&Wmnij, &Cme, &W, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);

    /** - C_n^b <Mn|Ij> **/
    dpd_buf4_init(&Wmnij, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");
    dpd_contract424(&Wmnij, &Cme, &W, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);

    /** - C_N^B <mN|iJ> **/
    dpd_buf4_init(&Wmnij, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WmBiJ (mB,iJ)");
    dpd_contract424(&Wmnij, &CME, &W, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);

    /* term 2: - P(ij) <mb||ei> Cje -> Wmbij */

    /** + P(IJ) C_J^E <MB||IE> -> WMBIJ **/
    dpd_buf4_init(&Z2, CC_TMP0, C_irr, 10, 0, 10, 0, 0, "Z1(MB,IJ)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_contract424(&C, &CME, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z2, CC_TMP1, pqsr, 10, 0, "Z2(MB,JI)");

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 10, 0, 10, 0, 0, "Z2(MB,JI)");
    dpd_buf4_axpy(&Z1, &Z2, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 0, 10, 2, 0, "HC1 WMBIJ (MB,I>J)");
    dpd_buf4_axpy(&Z2, &W, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&W);

    /** - P(ij) C_j^e ( <mb||ie> ) -> WMBIJ **/

    dpd_buf4_init(&Z2, CC_TMP0, C_irr, 10, 0, 10, 0, 0, "Z1(mb,ij)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_contract424(&C, &Cme, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z2, CC_TMP1, pqsr, 10, 0, "Z2(mb,ji)");

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 10, 0, 10, 0, 0, "Z2(mb,ji)");
    dpd_buf4_axpy(&Z1, &Z2, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 0, 10, 2, 0, "HC1 Wmbij (mb,i>j)");
    dpd_buf4_axpy(&Z2, &W, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&W);

    /* <Mb||Ie> Cje -> W(Mb,Ij) */
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_contract424(&C, &Cme, &W, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&C); 
                    
    /* CIE <Mb|Ej> -> W(Mb,Ij) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&CME, &D, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /** C_J^E <mB||iE> = Cje <mB|iE> **/
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WmBiJ (mB,iJ)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_contract424(&C, &CME, &W, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&C);

    /** -C_i^e <mB||Je> = +Cie <mB|eJ> = +Cie <mJ|eB> **/
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&Cme, &D, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* do purge before sort */
    purge_Wmbij(C_irr);

    /* do final sort to get (Ij,Mb) */

    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 2, 10, 2, 0, "HC1 WMBIJ (MB,I>J)");
    dpd_buf4_sort(&W, CC3_HC1, rspq, 2, 10, "HC1 WMBIJ (I>J,MB)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 2, 10, 2, 0, "HC1 Wmbij (mb,i>j)");
    dpd_buf4_sort(&W, CC3_HC1, rspq, 2, 10, "HC1 Wmbij (i>j,mb)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WMbIj (Mb,Ij)");
    dpd_buf4_sort(&W, CC3_HC1, rspq, 0, 10, "HC1 WMbIj (Ij,Mb)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HC1, C_irr, 10, 0, 10, 0, 0, "HC1 WmBiJ (mB,iJ)");
    dpd_buf4_sort(&W, CC3_HC1, rspq, 0, 10, "HC1 WmBiJ (iJ,mB)");
    dpd_buf4_close(&W);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /** - C_N^B W_MNIJ --> W(MB,IJ) **/
    dpd_buf4_init(&Wmnij, CC_AINTS, 0, 0, 2, 0, 0, 1, "A <IJ|KL>");
    dpd_buf4_init(&W, CC3_HC1, C_irr, 20, 2, 20, 2, 0, "HC1 WMBIJ (MB,I>J)");
    dpd_contract424(&Wmnij, &CME, &W, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);

    /** - C_n^b W_mnij **/
    dpd_buf4_init(&Wmnij, CC_AINTS, 0, 10, 12, 10, 10, 1, "A <ij|kl>");
    dpd_buf4_init(&W, CC3_HC1, C_irr, 30, 12, 30, 12, 0, "HC1 Wmbij (mb,i>j)");
    dpd_contract424(&Wmnij, &Cme, &W, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);

    /** - C_n^b W_MnIj **/
    dpd_buf4_init(&Wmnij, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    dpd_buf4_init(&W, CC3_HC1, C_irr, 24, 22, 24, 22, 0, "HC1 WMbIj (Mb,Ij)");
    dpd_contract424(&Wmnij, &Cme, &W, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);

    /** - C_N^B W_mNiJ **/
    dpd_buf4_init(&Wmnij, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    dpd_buf4_sort(&Wmnij, CC_TMP0, qpsr, 23, 23, "A <iJ|kL>");
    dpd_buf4_close(&Wmnij);
    dpd_buf4_init(&Wmnij, CC_TMP0, 0, 23, 23, 23, 23, 0, "A <iJ|kL>");
    dpd_buf4_init(&W, CC3_HC1, C_irr, 27, 23, 27, 23, 0, "HC1 WmBiJ (mB,iJ)");
    dpd_contract424(&Wmnij, &CME, &W, 1, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);

    /* term 2: - P(ij) <mb||ei> Cje -> Wmbij */

    /** + P(IJ) C_J^E <MB||IE> -> WMBIJ **/
    dpd_buf4_init(&Z2, CC_TMP0, C_irr, 20, 0, 20, 0, 0, "Z1(MB,IJ)");
    dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    dpd_contract424(&C, &CME, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z2, CC_TMP1, pqsr, 20, 0, "Z2(MB,JI)");

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 20, 0, 20, 0, 0, "Z2(MB,JI)");
    dpd_buf4_axpy(&Z1, &Z2, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&W, CC3_HC1, C_irr, 20, 0, 20, 2, 0, "HC1 WMBIJ (MB,I>J)");
    dpd_buf4_axpy(&Z2, &W, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&W);

    /** - P(ij) C_j^e ( <mb||ie> ) -> WMBIJ **/

    dpd_buf4_init(&Z2, CC_TMP0, C_irr, 30, 10, 30, 10, 0, "Z1(mb,ij)");
    dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_contract424(&C, &Cme, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z2, CC_TMP1, pqsr, 30, 10, "Z2(mb,ji)");

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 30, 10, 30, 10, 0, "Z2(mb,ji)");
    dpd_buf4_axpy(&Z1, &Z2, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&W, CC3_HC1, C_irr, 30, 10, 30, 12, 0, "HC1 Wmbij (mb,i>j)");
    dpd_buf4_axpy(&Z2, &W, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&W);

    /* <Mb||Ie> Cje -> W(Mb,Ij) */
    dpd_buf4_init(&W, CC3_HC1, C_irr, 24, 22, 24, 22, 0, "HC1 WMbIj (Mb,Ij)");
    dpd_buf4_init(&C, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    dpd_contract424(&C, &Cme, &W, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&C); 
                    
    /* CIE <Mb|Ej> -> W(Mb,Ij) */
    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_contract244(&CME, &D, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /** C_J^E <mB||iE> = C^J_E <mB|iE> **/
    dpd_buf4_init(&W, CC3_HC1, C_irr, 27, 23, 27, 23, 0, "HC1 WmBiJ (mB,iJ)");
    dpd_buf4_init(&C, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    dpd_contract424(&C, &CME, &W, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&C);

    /** -C_i^e <mB||Je> = +Cie <mB|eJ> = +Cie <mJ|eB> **/
    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_contract244(&Cme, &D, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /** do final sort to (Ij,Mb) **/
    dpd_buf4_init(&W, CC3_HC1, C_irr, 20, 2, 20, 2, 0, "HC1 WMBIJ (MB,I>J)");
    dpd_buf4_sort(&W, CC3_HC1, rspq, 2, 20, "HC1 WMBIJ (I>J,MB)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC3_HC1, C_irr, 30, 12, 30, 12, 0, "HC1 Wmbij (mb,i>j)");
    dpd_buf4_sort(&W, CC3_HC1, rspq, 12, 30, "HC1 Wmbij (i>j,mb)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC3_HC1, C_irr, 24, 22, 24, 22, 0, "HC1 WMbIj (Mb,Ij)");
    dpd_buf4_sort(&W, CC3_HC1, rspq, 22, 24, "HC1 WMbIj (Ij,Mb)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC3_HC1, C_irr, 27, 23, 27, 23, 0, "HC1 WmBiJ (mB,iJ)");
    dpd_buf4_sort(&W, CC3_HC1, rspq, 23, 27, "HC1 WmBiJ (iJ,mB)");
    dpd_buf4_close(&W);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
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

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);

    dpd_buf4_init(&WMbEj, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &CME, &WMbEj, 3, 1, 0, 1, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&WMbEj);

    
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_contract424(&F, &CME, &Z, 3, 1, 0, -1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_sort(&Z, CC_TMP0, qpsr, 10, 10, "WMbeJ"); /* (Mb,Je) */
    dpd_buf4_close(&Z);

    /* - Cnb <mn||ej> -> Wmbej */

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&WMbEj, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&E, &CME, &WMbEj, 3, 0, 1, -1, 1.0);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&WMbeJ, CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ"); 
    dpd_contract424(&E, &CME, &WMbeJ, 1, 0, 1, 1, 1.0);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_close(&E);

    dpd_file2_close(&CME);

    /* Sort to (ME,JB) */

    dpd_buf4_init(&WMbEj, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    dpd_buf4_sort(&WMbEj, CC3_HC1, prsq, 10, 10, "HC1 WMbEj (ME,jb)");
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WMbeJ, CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_sort(&WMbeJ, CC3_HC1, psrq, 10, 10, "HC1 WMbeJ (Me,Jb)");
    dpd_buf4_close(&WMbeJ);

  }
  else if(params.ref == 1) { /** ROHF **/


    /* + Cjf <mb||ef> -> Wmbej*/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&WMBEJ, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMBEJ");
    dpd_contract424(&F, &CME, &WMBEJ, 3, 1, 0, 1, 0);
    dpd_buf4_close(&WMBEJ);
    dpd_buf4_init(&Wmbej, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "Wmbej");
    dpd_contract424(&F, &Cme, &Wmbej, 3, 1, 0, 1, 0);
    dpd_buf4_close(&Wmbej);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&WMbEj, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&F, &Cme, &WMbEj, 3, 1, 0, 1, 0);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_init(&WmBeJ, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WmBeJ");
    dpd_contract424(&F, &CME, &WmBeJ, 3, 1, 0, 1, 0);
    dpd_buf4_close(&WmBeJ);

    dpd_buf4_init(&WMbeJ, CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract244(&CME, &F, &WMbeJ, 1, 2, 1, -1, 0);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_init(&WmBEj, CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WmBEj");
    dpd_contract244(&Cme, &F, &WmBEj, 1, 2, 1, -1, 0);
    dpd_buf4_close(&WmBEj);
    dpd_buf4_close(&F);

    /* - Cnb <mn||ej> -> Wmbej */

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    dpd_buf4_init(&WMBEJ, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMBEJ");
    dpd_contract424(&E, &CME, &WMBEJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMBEJ);
    dpd_buf4_init(&Wmbej, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "Wmbej");
    dpd_contract424(&E, &Cme, &Wmbej, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Wmbej);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&WMbEj, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    dpd_contract424(&E, &Cme, &WMbEj, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_init(&WmBeJ, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WmBeJ");
    dpd_contract424(&E, &CME, &WmBeJ, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WmBeJ);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&WMbeJ, CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract424(&E, &Cme, &WMbeJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_init(&WmBEj, CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WmBEj");
    dpd_contract424(&E, &CME, &WmBEj, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WmBEj);
    dpd_buf4_close(&E);

    /* Convert to (ME,JB) for remaining terms */

    dpd_buf4_init(&WMBEJ, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMBEJ");
    dpd_buf4_sort(&WMBEJ, CC3_HC1, prsq, 10, 10, "HC1 WMBEJ (ME,JB)");
    dpd_buf4_close(&WMBEJ);

    dpd_buf4_init(&Wmbej, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "Wmbej");
    dpd_buf4_sort(&Wmbej, CC3_HC1, prsq, 10, 10, "HC1 Wmbej (me,jb)");
    dpd_buf4_close(&Wmbej);

    dpd_buf4_init(&WMbEj, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WMbEj");
    dpd_buf4_sort(&WMbEj, CC3_HC1, prsq, 10, 10, "HC1 WMbEj (ME,jb)");
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WmBeJ, CC_TMP0, C_irr, 10, 11, 10, 11, 0, "WmBeJ");
    dpd_buf4_sort(&WmBeJ, CC3_HC1, prsq, 10, 10, "HC1 WmBeJ (me,JB)");
    dpd_buf4_close(&WmBeJ);

    dpd_buf4_init(&WMbeJ, CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_sort(&WMbeJ, CC3_HC1, psrq, 10, 10, "HC1 WMbeJ (Me,Jb)");
    dpd_buf4_close(&WMbeJ);

    dpd_buf4_init(&WmBEj, CC_TMP0, C_irr, 10, 10, 10, 10, 0, "WmBEj");
    dpd_buf4_sort(&WmBEj, CC3_HC1, psrq, 10, 10, "HC1 WmBEj (mE,jB)");
    dpd_buf4_close(&WmBEj);

  } /** ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* F -> Wmbej */ /* + Cjf <mb||ef> -> Wmbej*/

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 20, 21, 20, 21, 0, "WMBEJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_contract424(&F, &CME, &W, 3, 1, 0, 1, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 30, 31, 30, 31, 0, "Wmbej");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_contract424(&F, &Cme, &W, 3, 1, 0, 1, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 24, 26, 24, 26, 0, "WMbEj");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract424(&F, &Cme, &W, 3, 1, 0, 1, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 27, 25, 27, 25, 0, "WmBeJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract424(&F, &CME, &W, 3, 1, 0, 1, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract244(&CME, &F, &W, 1, 2, 1, -1, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract244(&Cme, &F, &W, 1, 2, 1, -1, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    /* - Cnb <mn||ej> -> Wmbej */

    dpd_buf4_init(&W, CC_TMP0, C_irr, 20, 21, 20, 21, 0, "WMBEJ");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    dpd_contract424(&E, &CME, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 30, 31, 30, 31, 0, "Wmbej");
    dpd_buf4_init(&E, CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    dpd_contract424(&E, &Cme, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 24, 26, 24, 26, 0, "WMbEj");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_contract424(&E, &Cme, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 27, 25, 27, 25, 0, "WmBeJ");
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    dpd_contract424(&E, &CME, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_contract424(&E, &Cme, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_contract424(&E, &CME, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    /* Convert to (ME,JB) for remaining terms */

    dpd_buf4_init(&W, CC_TMP0, C_irr, 20, 21, 20, 21, 0, "WMBEJ");
    dpd_buf4_sort(&W, CC3_HC1, prsq, 20, 20, "HC1 WMBEJ (ME,JB)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 30, 31, 30, 31, 0, "Wmbej");
    dpd_buf4_sort(&W, CC3_HC1, prsq, 30, 30, "HC1 Wmbej (me,jb)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 24, 26, 24, 26, 0, "WMbEj");
    dpd_buf4_sort(&W, CC3_HC1, prsq, 20, 30, "HC1 WMbEj (ME,jb)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 27, 25, 27, 25, 0, "WmBeJ");
    dpd_buf4_sort(&W, CC3_HC1, prsq, 30, 20, "HC1 WmBeJ (me,JB)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_buf4_sort(&W, CC3_HC1, psrq, 24, 24, "HC1 WMbeJ (Me,Jb)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 27, 27, 27, 27, 0, "WmBEj");
    dpd_buf4_sort(&W, CC3_HC1, psrq, 27, 27, "HC1 WmBEj (mE,jB)");
    dpd_buf4_close(&W);
  }

  return;
}



void HC1_Wabei(int i, int C_irr) { 
  dpdbuf4 Z, Z1, Z2, Z3;
  dpdbuf4 B, C, D, E, F, W;
  dpdfile2 Cme, CME;
  char CME_lbl[32], Cme_lbl[32];

  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "CC3 Z(Ab,Ei)");
    dpd_buf4_init(&Z2, CC_TMP0, C_irr, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");

    /* Z1(Ab,Ei) <-- <Ab|Ef> * C(i,f) */
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract424(&B, &CME, &Z1, 3, 1, 0, 1, 0);
    dpd_buf4_close(&B);

    /* Z1(Ab,Ei) <--  - C(M,A) * <Mb|Ei> */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&CME, &D, &Z1, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);

    /* Z2(Ei,Ab) <-- - <mA,iE> C(m,b) */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_sort(&C, CC_TMP0, qpsr, 11, 11, "CC3 Z(Ei,Am)");
    dpd_buf4_close(&C);

    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "CC3 Z(Ei,Am)");
    dpd_contract424(&Z, &CME, &Z2, 3, 0, 0, -1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_close(&Z2);

    /* W(Ab,Ei) = Z1(Ab,Ei) + Z2(Ei,Ab) */
    dpd_buf4_sort_axpy(&Z1, CC_TMP0, rspq, 11, 5, "CC3 Z(Ei,Ab)", 1);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z2, CC_TMP0, C_irr, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");
    dpd_buf4_sort(&Z2, CC3_HC1, qpsr, 10, 5, "CC3 WAbEi (Ie,Ab)");

    dpd_file2_close(&CME);
  }

  else if (params.ref == 1) { /* ROHF */
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 7, 11, 7, 11, 0, "WABEI (A>B,EI)");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 5, 5, 5, 1, "B <ab|cd>");
    dpd_contract424(&B, &CME, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 7, 11, 7, 11, 0, "Wabei (a>b,ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 5, 5, 5, 1, "B <ab|cd>");
    dpd_contract424(&B, &Cme, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WAbEi (Ab,Ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract424(&B, &Cme, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WaBeI (aB,eI)");
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract424(&B, &CME, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z1);

    /** -CMA <MI||EB> + CMB <MI||EA> **/

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 5, 11, 5, 11, 0, "Z (AB,EI)");

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    dpd_contract244(&CME, &C, &Z1, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z1, CC_TMP1, qprs, 5, 11, "Z (BA,EI)");

    dpd_buf4_init(&W, CC_TMP0, C_irr, 5, 11, 7, 11, 0, "WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z1, &W, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 5, 11, 5, 11, 0, "Z (BA,EI)");
    dpd_buf4_axpy(&Z1, &W, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&W);

    /** -Cma <mi||eb> + Cmb <mi||ea> **/

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 5, 11, 5, 11, 0, "Z (ab,ei)");

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    dpd_contract244(&Cme, &C, &Z1, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z1, CC_TMP1, qprs, 5, 11, "Z (ba,ei)");

    dpd_buf4_init(&W, CC_TMP0, C_irr, 5, 11, 7, 11, 0, "Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z1, &W, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 5, 11, 5, 11, 0, "Z (ba,ei)");
    dpd_buf4_axpy(&Z1, &W, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&W);

    /** -CMA <Mi||Eb> - Cmb <mA||iE> **/

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WAbEi (Ab,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&CME, &D, &Z1, 0, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 5, 11, 5, 11, 0, "WAbEi (bA,Ei)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia|jb> (ia,bj)");
    dpd_contract244(&Cme, &C, &Z1, 0, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 5, 11, "WAbEi (Ab,Ei)", 1.0);
    dpd_buf4_close(&Z1);

    /** **/
    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WaBeI (aB,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&Cme, &D, &Z1, 0, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 5, 11, 5, 11, 0, "WaBeI (Ba,eI)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia|jb> (ia,bj)");
    dpd_contract244(&CME, &C, &Z1, 0, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 5, 11, "WaBeI (aB,eI)", 1.0);
    dpd_buf4_close(&Z1);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);

    /* final sort to (EI,AB) */

    dpd_buf4_init(&W, CC_TMP0, C_irr, 5, 11, 7, 11, 0, "WABEI (A>B,EI)");
    dpd_buf4_sort(&W, CC_TMP0, rspq, 11, 7, "WABEI (EI,A>B)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, C_irr, 5, 11, 7, 11, 0, "Wabei (a>b,ei)");
    dpd_buf4_sort(&W, CC_TMP0, rspq, 11, 7, "Wabei (ei,a>b)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WAbEi (Ab,Ei)");
    dpd_buf4_sort(&W, CC_TMP0, rspq, 11, 5, "WAbEi (Ei,Ab)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, C_irr, 5, 11, 5, 11, 0, "WaBeI (aB,eI)");
    dpd_buf4_sort(&W, CC_TMP0, rspq, 11, 5, "WaBeI (eI,aB)");
    dpd_buf4_close(&W);

    purge_Wabei(C_irr);

    /* final sort to Wabei (IE,AB) */
    dpd_buf4_init(&W, CC_TMP0, C_irr, 11, 7, 11, 7, 0, "WABEI (EI,A>B)");
    dpd_buf4_sort(&W, CC3_HC1, qprs, 10, 7, "HC1 WABEI (IE,A>B)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, C_irr, 11, 7, 11, 7, 0, "Wabei (ei,a>b)");
    dpd_buf4_sort(&W, CC3_HC1, qprs, 10, 7, "HC1 Wabei (ie,a>b)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, C_irr, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
    dpd_buf4_sort(&W, CC3_HC1, qprs, 10, 5, "HC1 WAbEi (iE,Ab)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, C_irr, 11, 5, 11, 5, 0, "WaBeI (eI,aB)");
    dpd_buf4_sort(&W, CC3_HC1, qprs, 10, 5, "HC1 WaBeI (Ie,aB)");
    dpd_buf4_close(&W);
  }
  else if (params.ref == 2) { /* UHF */
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    /* term 1, Cif <ab||ef> */

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 7, 21, 7, 21, 0, "WABEI (A>B,EI)");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
    dpd_contract424(&B, &CME, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 17, 31, 17, 31, 0, "Wabei (a>b,ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 17, 15, 15, 15, 1, "B <ab|cd>");
    dpd_contract424(&B, &Cme, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 28, 26, 28, 26, 0, "WAbEi (Ab,Ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_contract424(&B, &Cme, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z, CC_TMP0, C_irr, 24, 28, 24, 28, 0, "Z(Ie,Ba)");
    dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_contract244(&CME, &B, &Z, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&B);
    /** Z(Ie,Ba) --> W'(aB,eI) **/
    /* srqp seems to have a bug
     dpd_buf4_sort(&Z, CC_TMP0, srqp, 29, 25, "WaBeI (aB,eI)");
    */
    dpd_buf4_sort(&Z, CC_TMP0, rspq, 28, 24, "WaBeI (Ba,Ie) 1");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, C_irr, 28, 24, 28, 24, 0, "WaBeI (Ba,Ie)");
    dpd_buf4_sort(&Z, CC_TMP0, qpsr, 29, 25, "WaBeI (aB,eI)");
    dpd_buf4_close(&Z);

    /** UHF term 2 **/
    /* -CMA <MI||EB> + CMB <MI||EA> **/

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 5, 21, 5, 21, 0, "Z (AB,EI)");

    dpd_buf4_init(&C, CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    dpd_contract244(&CME, &C, &Z1, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z1, CC_TMP1, qprs, 5, 21, "Z (BA,EI)");

    dpd_buf4_init(&W, CC_TMP0, C_irr, 5, 21, 7, 21, 0, "WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z1, &W, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 5, 21, 5, 21, 0, "Z (BA,EI)");
    dpd_buf4_axpy(&Z1, &W, -1.0);
    dpd_buf4_close(&Z1);

    /** -Cma <mi||eb> + Cmb <mi||ea> **/

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 15, 31, 15, 31, 0, "Z (ab,ei)");

    dpd_buf4_init(&C, CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    dpd_contract244(&Cme, &C, &Z1, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z1, CC_TMP1, qprs, 15, 31, "Z (ba,ei)");

    dpd_buf4_init(&W, CC_TMP0, C_irr, 15, 31, 17, 31, 0, "Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z1, &W, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 15, 31, 15, 31, 0, "Z (ba,ei)");
    dpd_buf4_axpy(&Z1, &W, -1.0);
    dpd_buf4_close(&Z1);

    /** -CMA <Mi||Eb> - Cmb <mA|iE> **/

    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 28, 26, 28, 26, 0, "WAbEi (Ab,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_contract244(&CME, &D, &Z1, 0, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 29, 26, 29, 26, 0, "WAbEi (bA,Ei)");
    dpd_buf4_init(&C, CC_CINTS, 0, 27, 26, 27, 26, 0, "C <Ai|Bj> (iA,Bj)");
    dpd_contract244(&Cme, &C, &Z1, 0, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 28, 26, "WAbEi (Ab,Ei)", 1.0);
    dpd_buf4_close(&Z1);

    /**  **/
    dpd_buf4_init(&Z1, CC_TMP0, C_irr, 29, 25, 29, 25, 0, "WaBeI (aB,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_contract244(&Cme, &D, &Z1, 0, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z1, CC_TMP1, C_irr, 28, 25, 28, 25, 0, "WaBeI (Ba,eI)");
    dpd_buf4_init(&C, CC_CINTS, 0, 24, 25, 24, 25, 0, "C <Ia|Jb> (Ia,bJ)");
    dpd_contract244(&CME, &C, &Z1, 0, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&C);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 29, 25, "WaBeI (aB,eI)", 1.0);
    dpd_buf4_close(&Z1);

    /* final sort and storage */
    dpd_buf4_init(&W, CC_TMP0, C_irr, 7, 21, 7, 21, 0, "WABEI (A>B,EI)");
    dpd_buf4_sort(&W, CC_TMP0, rspq, 21, 7, "WABEI (EI,A>B)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 21, 7, 21, 7, 0, "WABEI (EI,A>B)");
    dpd_buf4_sort(&W, CC3_HC1, qprs, 20, 7, "HC1 WABEI (IE,B>A)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC3_HC1, C_irr, 20, 7, 20, 7, 0, "HC1 WABEI (IE,B>A)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    /* final sort and storage */
    dpd_buf4_init(&W, CC_TMP0, C_irr, 17, 31, 17, 31, 0, "Wabei (a>b,ei)");
    dpd_buf4_sort(&W, CC_TMP0, rspq, 31, 17, "Wabei (ei,a>b)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, C_irr, 31, 17, 31, 17, 0, "Wabei (ei,a>b)");
    dpd_buf4_sort(&W, CC3_HC1, qprs, 30, 17, "HC1 Wabei (ie,b>a)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC3_HC1, C_irr, 30, 17, 30, 17, 0, "HC1 Wabei (ie,b>a)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    /* final sort and storage */
    dpd_buf4_init(&W, CC_TMP0, C_irr, 28, 26, 28, 26, 0, "WAbEi (Ab,Ei)");
    dpd_buf4_sort(&W, CC_TMP0, rspq, 26, 28, "WAbEi (Ei,Ab)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, C_irr, 26, 28, 26, 28, 0, "WAbEi (Ei,Ab)");
    dpd_buf4_sort(&W, CC3_HC1, qpsr, 27, 29, "HC1 WAbEi (iE,bA)");
    dpd_buf4_close(&W);

    /* final sort and storage */
    dpd_buf4_init(&W, CC_TMP0, C_irr, 29, 25, 29, 25, 0, "WaBeI (aB,eI)");
    dpd_buf4_sort(&W, CC_TMP0, rspq, 25, 29, "WaBeI (eI,aB)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, C_irr, 25, 29, 25, 29, 0, "WaBeI (eI,aB)");
    dpd_buf4_sort(&W, CC3_HC1, qpsr, 24, 28, "HC1 WaBeI (Ie,Ba)");
    dpd_buf4_close(&W);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
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
  dpd_file2_init(&FME, CC3_HC1, C_irr, 0, 1, "HC1 FME");
  dpd_file2_mat_init(&FME);
  dpd_file2_mat_rd(&FME);
  for(h=0; h < nirreps; h++) {
    for(m=0; m<occpi[h]; m++)
      for(e=(virtpi[h]-openpi[h]); e<virtpi[h]; e++)
	FME.matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(&FME);
  dpd_file2_mat_close(&FME);
  dpd_file2_close(&FME);

  /* Purge Fme matrix elements */
  dpd_file2_init(&Fme, CC3_HC1, C_irr, 0, 1, "HC1 Fme");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);
  for(h=0; h < nirreps; h++) {
    for(e=0; e<virtpi[h]; e++)
      for(m=(occpi[h]-openpi[h]); m<occpi[h]; m++)
	Fme.matrix[h][m][e] = 0.0;
  }
  dpd_file2_mat_wrt(&Fme);
  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  /* Purge Wmnij matrix elements */
  dpd_file4_init(&W, CC3_HC1, C_irr, 2, 2,"HC1 Wmnij (m>n,i>j)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HC1, C_irr, 0, 0,"HC1 WMnIj (Mn,Ij)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  /* Purge Wmbej matrix elements */
  dpd_file4_init(&W, CC3_HC1, C_irr, 10, 10,"HC1 WMBEJ (ME,JB)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  dpd_file4_init(&W, CC3_HC1, C_irr, 10, 10,"HC1 Wmbej (me,jb)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  dpd_file4_init(&W, CC3_HC1, C_irr, 10, 10,"HC1 WMbEj (ME,jb)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  dpd_file4_init(&W, CC3_HC1, C_irr, 10, 10,"HC1 WmBeJ (me,JB)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  dpd_file4_init(&W, CC3_HC1, C_irr, 10, 10,"HC1 WmBEj (mE,jB)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);


  /* Purge Wamef matrix elements */
  dpd_file4_init(&W, CC3_HC1, C_irr, 11, 7,"HC1 WAMEF (AM,E>F)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HC1, C_irr, 11, 7,"HC1 Wamef (am,e>f)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      m = W.params->roworb[h][ma][1];
      msym = W.params->qsym[m];
      M = m - occ_off[msym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	if (M >=  (occpi[msym] - openpi[msym]))
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HC1, C_irr, 11, 5,"HC1 WAmEf (Am,Ef)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HC1, C_irr, 11, 5,"HC1 WaMeF (aM,eF)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      for(ef=0; ef< W.params->coltot[h]; ef++) {
	f = W.params->colorb[h][ef][1];
	fsym = W.params->ssym[f];
	F = f - vir_off[fsym];
	if (F >= (virtpi[fsym] - openpi[fsym]))
	  W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);



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

  dpd_file4_init(&W, CC3_HC1, C_irr, 0, 11,"HC1 WMnIe (Mn,eI)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      n = W.params->roworb[h][mn][1];
      nsym = W.params->qsym[n];
      N = n - occ_off[nsym];
      for(ei=0; ei<W.params->coltot[h]; ei++) {
	if (N >= (occpi[nsym] - openpi[nsym]))
	  W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }

  dpd_file4_init(&W, CC3_HC1, C_irr, 2, 11, "HC1 WMNIE (M>N,EI)");
  for(h=0; h < W.params->nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mn=0; mn<W.params->rowtot[h]; mn++) {
      for(ei=0; ei<W.params->coltot[h]; ei++) {
        e = W.params->colorb[h][ei][0];
        esym = W.params->rsym[e];
        E = e - vir_off[esym];
        if (E >= (virtpi[esym] - openpi[esym]))
          W.matrix[h][mn][ei] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HC1, C_irr, 2, 11,"HC1 Wmnie (m>n,ei)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HC1, C_irr, 0, 11,"HC1 WmNiE (mN,Ei)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);
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

  dpd_file4_init(&W, CC3_HC1, C_irr, 10, 2,"HC1 WMBIJ (MB,I>J)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      b = W.params->roworb[h][mb][1];
      bsym = W.params->qsym[b];
      B = b - vir_off[bsym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	if (B >= (virtpi[bsym] - openpi[bsym]))
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HC1, C_irr, 10, 2,"HC1 Wmbij (mb,i>j)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HC1, C_irr, 10, 0,"HC1 WMbIj (Mb,Ij)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      for(ij=0; ij<W.params->coltot[h]; ij++) {
	j = W.params->colorb[h][ij][1];
	jsym = W.params->ssym[j];
	J = j - occ_off[jsym];
	if (J >= (occpi[jsym] - openpi[jsym]))
	  W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HC1, C_irr, 10, 0,"HC1 WmBiJ (mB,iJ)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);
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

  dpd_file4_init(&W, CC_TMP0, C_irr, 11, 7,"WABEI (EI,A>B)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_TMP0, C_irr, 11, 7,"Wabei (ei,a>b)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      i = W.params->roworb[h][ei][1];
      isym = W.params->qsym[i];
      I = i - occ_off[isym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	if (I >= (occpi[isym] - openpi[isym]))
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_TMP0, C_irr, 11, 5,"WAbEi (Ei,Ab)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
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
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_TMP0, C_irr, 11, 5,"WaBeI (eI,aB)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      for(ab=0; ab<W.params->coltot[h]; ab++) {
	b = W.params->colorb[h][ab][1];
	bsym = W.params->ssym[b];
	B = b - vir_off[bsym];
	if (B >= (virtpi[bsym] - openpi[bsym]))
	  W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);
}

}} // namespace psi::cceom
