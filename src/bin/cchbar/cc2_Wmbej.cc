/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* Wmbej_build(): Computes all contributions to the Wmbej HBAR matrix
** elements.  These are defined in terms of spin orbitals as:
**
** Wmbej = <mb||ej> + t_j^f <mb||ef> - t_n^b <mn||ej> 
**         - { t_jn^fb + t_j^f t_n^b } <mn||ef>
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
**
** There are Wmbej six spin cases, which are stored and named
** as follows:
**
** Spin Case    Storage    Name
** ----------   ---------  -------
** WMBEJ        (ME,JB)    "WMBEJ"
** Wmbej        (me,jb)    "Wmbej"
** WMbEj        (ME,jb)    "WMbEj"
** WmBeJ        (me,JB)    "WmBeJ"
** WMbeJ        (Me,bJ)    "WMbeJ"
** WmBEj        (mE,Bj)    "WmBEj"
** -------------------------------
**
** TDC, June 2002
*/

void cc2_Wmbej_build(void) {
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ;
  dpdbuf4 tIAJB, tjAIb, tiajb, tIAjb, tiaJB, tIbjA;
  dpdbuf4 D, C, F, E, X, Y, Y1, t2, W, Z;
  dpdfile2 tIA, tia;

  if(params.ref == 0) { /** RHF **/

    /* <mb||ej> -> Wmbej */

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 WMbeJ (Mb,Je)", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "CC2 WMbEj");
    dpd_buf4_close(&D);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    /* F -> Wmbej */
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    dpd_contract424(&F, &tIA, &WMbEj, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, -1, 0);
    dpd_buf4_sort(&Z, CC_TMP0, qpsr, 10, 10, "Z(Mb,Je)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Je)");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Mb,Je)");
    dpd_buf4_axpy(&Z, &WMbeJ, 1.0);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&F);

    /* E -> Wmbej */
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    dpd_contract424(&E, &tIA, &WMbEj, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Mb,Je)");
    dpd_contract424(&E, &tIA, &WMbeJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA);


    /* Sort to (ME,JB) */
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    dpd_buf4_sort(&WMbEj, CC2_HET1, prsq, 10, 10, "CC2 WMbEj (ME,jb)");
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Mb,Je)");
    dpd_buf4_sort(&WMbeJ, CC2_HET1, psrq, 10, 10, "CC2 WMbeJ (Me,Jb)");
    dpd_buf4_close(&WMbeJ);


    /* T1^2 -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    /*** ABAB ***/
    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj (ME,jb)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABBA ***/
    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Me,Jb)");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);

  }

  else if(params.ref == 1) { /** ROHF **/

    /* W(mb,je) <-- <mb||ej> */

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 WMBEJ", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 Wmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 WmBEj", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 WMbeJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "CC2 WMbEj");
    dpd_buf4_copy(&D, CC_TMP0, "CC2 WmBeJ");
    dpd_buf4_close(&D);

    /* F -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMBEJ");
    dpd_contract424(&F, &tIA, &WMBEJ, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WMBEJ);
    dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Wmbej");
    dpd_contract424(&F, &tia, &Wmbej, 3, 1, 0, 1, 1);
    dpd_buf4_close(&Wmbej);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    dpd_contract424(&F, &tia, &WMbEj, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WmBeJ");
    dpd_contract424(&F, &tIA, &WmBeJ, 3, 1, 0, 1, 1);
    dpd_buf4_close(&WmBeJ);

    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    dpd_contract244(&tIA, &F, &WMbeJ, 1, 2, 1, -1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
    dpd_contract244(&tia, &F, &WmBEj, 1, 2, 1, -1, 1);
    dpd_buf4_close(&WmBEj);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /* E -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMBEJ");
    dpd_contract424(&E, &tIA, &WMBEJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMBEJ);
    dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Wmbej");
    dpd_contract424(&E, &tia, &Wmbej, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Wmbej);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    dpd_contract424(&E, &tia, &WMbEj, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WMbEj);
    dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WmBeJ");
    dpd_contract424(&E, &tIA, &WmBeJ, 3, 0, 1, -1, 1);
    dpd_buf4_close(&WmBeJ);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    dpd_contract424(&E, &tia, &WMbeJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WMbeJ);
    dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
    dpd_contract424(&E, &tIA, &WmBEj, 1, 0, 1, 1, 1);
    dpd_buf4_close(&WmBEj);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA);  
    dpd_file2_close(&tia);

    /* Convert to (ME,JB) for remaining terms */

    dpd_buf4_init(&WMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMBEJ");
    dpd_buf4_sort(&WMBEJ, CC2_HET1, prsq, 10, 10, "CC2 WMBEJ");
    dpd_buf4_close(&WMBEJ);

    dpd_buf4_init(&Wmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Wmbej");
    dpd_buf4_sort(&Wmbej, CC2_HET1, prsq, 10, 10, "CC2 Wmbej");
    dpd_buf4_close(&Wmbej);

    dpd_buf4_init(&WMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    dpd_buf4_sort(&WMbEj, CC2_HET1, prsq, 10, 10, "CC2 WMbEj");
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&WmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WmBeJ");
    dpd_buf4_sort(&WmBeJ, CC2_HET1, prsq, 10, 10, "CC2 WmBeJ");
    dpd_buf4_close(&WmBeJ);

    dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    dpd_buf4_sort(&WMbeJ, CC2_HET1, psrq, 10, 10, "CC2 WMbeJ");
    dpd_buf4_close(&WMbeJ);

    dpd_buf4_init(&WmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
    dpd_buf4_sort(&WmBEj, CC2_HET1, psrq, 10, 10, "CC2 WmBEj");
    dpd_buf4_close(&WmBEj);

    /* X -> Wmbej */
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /*** AAAA ***/

    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMBEJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMBEJ");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BBBB ***/

    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 Wmbej");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 Wmbej");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABAB ***/
  
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1); 
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1); 
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BABA ***/

    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBeJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBeJ");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABBA ***/
  
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BAAB ***/

    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&t2);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);


    /* D(me,nf) * T2(jb,nf) --> W(me,jb) */
    dpd_buf4_init(&Wmbej, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 Wmbej");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tiajb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &tiajb, &Wmbej, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiajb);

    /* D(me,NF) * T2(jb,NF) --> W(me,jb) */
    dpd_buf4_init(&tiaJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_contract444(&D, &tiaJB, &Wmbej, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiaJB);

    dpd_buf4_close(&Wmbej);


    /* D(ME,NF) * T2(JB,NF) --> W(ME,JB) */
    dpd_buf4_init(&WMBEJ, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMBEJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tIAJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &tIAJB, &WMBEJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAJB);

    /* D(ME,nf) * T2(JB,nf) --> W(ME,JB) */
    dpd_buf4_init(&tIAjb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_contract444(&D, &tIAjb, &WMBEJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAjb);

    dpd_buf4_close(&WMBEJ);


    /* D(me,nf) * T2(JB,nf) --> W(me,JB) */
    dpd_buf4_init(&WmBeJ, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBeJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tIAjb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &tIAjb, &WmBeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAjb);

    /* D(me,NF) * T2(JB,NF) --> W(me,JB) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&tIAJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &tIAJB, &WmBeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAJB);

    dpd_buf4_close(&WmBeJ);


    /* D(ME,NF) * T2(jb,NF) --> W(ME,jb) */
    dpd_buf4_init(&WMbEj, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tiaJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&D, &tiaJB, &WMbEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiaJB);

    /* D(ME,nf) * T2(jb,nf) --> W(ME,jb) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&tiajb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &tiajb, &WMbEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiajb);

    dpd_buf4_close(&WMbEj);


    /* D(mE,Nf) * T2(jB,Nf) --> W(mE,jB) */
    dpd_buf4_init(&WmBEj, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&tjAIb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_contract444(&D, &tjAIb, &WmBEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tjAIb);

    dpd_buf4_close(&WmBEj);


    /* D(Me,nF) * T2(Jb,nF) --> W(Me,Jb) */
    dpd_buf4_init(&WMbeJ, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&tIbjA, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_contract444(&D, &tIbjA, &WMbeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIbjA);

    dpd_buf4_close(&WMbeJ);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* W(mb,je) <-- <mb||ej> */

    dpd_buf4_init(&C, CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 WMBEJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 Wmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_buf4_scmcopy(&D, CC_TMP0, "CC2 WMbEj", -1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_buf4_scmcopy(&D, CC_TMP0, "CC2 WmBeJ", -1);
    dpd_buf4_close(&D);

    /* F -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&W, CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 WMBEJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Wmbej");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 WMbEj");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract424(&F, &tia, &W, 3, 1, 0, -1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 WmBeJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract424(&F, &tIA, &W, 3, 1, 0, -1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /* D & E -> Wmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /*** AAAA ***/

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    dpd_buf4_scmcopy(&E, CC_TMP0, "YMNEJ", -1);
    dpd_buf4_close(&E);

    dpd_buf4_init(&Y, CC_TMP0, 0, 0, 21, 0, 21, 0, "YMNEJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_contract424(&D, &tIA, &Y, 3, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 WMBEJ");
    dpd_contract424(&Y, &tIA, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&Y);
    dpd_buf4_close(&W);

    /*** BBBB ***/

    dpd_buf4_init(&E, CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    dpd_buf4_scmcopy(&E, CC_TMP0, "Ymnej", -1);
    dpd_buf4_close(&E);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 31, 10, 31, 0, "Ymnej");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D, &tia, &Y, 3, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Wmbej");
    dpd_contract424(&Y, &tia, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&Y);
    dpd_buf4_close(&W);

    /*** ABAB ***/

    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_buf4_copy(&E, CC_TMP0, "YMnEj");
    dpd_buf4_close(&E);

    dpd_buf4_init(&Y, CC_TMP0, 0, 22, 26, 22, 26, 0, "YMnEj");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &tia, &Y, 3, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 WMbEj");
    dpd_contract424(&Y, &tia, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Y);
    dpd_buf4_close(&W);

    /*** BABA ***/

    dpd_buf4_init(&E, CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    dpd_buf4_copy(&E, CC_TMP0, "YmNeJ");
    dpd_buf4_close(&E);

    dpd_buf4_init(&Y, CC_TMP0, 0, 23, 25, 23, 25, 0, "YmNeJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract424(&D, &tIA, &Y, 3, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 WmBeJ");
    dpd_contract424(&Y, &tIA, &W, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Y);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /** Sort **/
    dpd_buf4_init(&W, CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 WMBEJ");
    dpd_buf4_sort(&W, CC2_HET1, prqs, 20, 21, "CC2 WMBEJ (ME,BJ)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Wmbej");
    dpd_buf4_sort(&W, CC2_HET1, prqs, 30, 31, "CC2 Wmbej (me,bj)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 WMbEj");
    dpd_buf4_sort(&W, CC2_HET1, prqs, 20, 31, "CC2 WMbEj (ME,bj)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 WmBeJ");
    dpd_buf4_sort(&W, CC2_HET1, prqs, 30, 21, "CC2 WmBeJ (me,BJ)");
    dpd_buf4_close(&W);

  } /** UHF **/

  return;
}

}} // namespace psi::cchbar
