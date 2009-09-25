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

/* Zmbej_build(): Computes all contributions to the Zmbej HBAR matrix
** elements.  These are defined in terms of spin orbitals as:
**
** Zmbej = <mb||ej> + t_j^f <mb||ef> - t_n^b <mn||ej> 
**         - { t_jn^fb + t_j^f t_n^b } <mn||ef>
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
**
** There are Zmbej six spin cases, which are stored and named
** as follows:
**
** Spin Case    Storage    Name
** ----------   ---------  -------
** ZMBEJ        (ME,JB)    "ZMBEJ"
** Zmbej        (me,jb)    "Zmbej"
** ZMbEj        (ME,jb)    "ZMbEj"
** ZmBeJ        (me,JB)    "ZmBeJ"
** ZMbeJ        (Me,bJ)    "ZMbeJ"
** ZmBEj        (mE,Bj)    "ZmBEj"
** -------------------------------
**
** TDC, June 2002
*/

void cc2_Zmbej_build(void) {
  dpdbuf4 ZMBEJ, Zmbej, ZMbEj, ZmBeJ, ZmBEj, ZMbeJ;
  dpdbuf4 tIAJB, tjAIb, tiajb, tIAjb, tiaJB, tIbjA;
  dpdbuf4 D, C, F, E, X, Y, Y1, t2, W, Z;
  dpdfile2 tIA, tia;

  if(params.ref == 0) { /** RHF **/

    /* <mb||ej> -> Zmbej */

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "CC2 ZMbEj");
    dpd_buf4_close(&D);

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 ZMbeJ (Mb,Je)", 1);
    dpd_buf4_close(&C);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    /* F -> Zmbej */
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&ZMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    dpd_contract424(&F, &tIA, &ZMbEj, 3, 1, 0, 1, 1);
    dpd_buf4_close(&ZMbEj);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, -1, 0);
    dpd_buf4_sort(&Z, CC_TMP0, qpsr, 10, 10, "Z(Mb,Je)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Je)");
    dpd_buf4_init(&ZMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ (Mb,Je)");
    dpd_buf4_axpy(&Z, &ZMbeJ, -1.0);
    dpd_buf4_close(&ZMbeJ);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&F);

    /* E -> Zmbej */

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&ZMbEj, CC_TMP0, 0, 11, 10, 11, 10, 0, "CC2 ZMbEj (Ej,Mb)");
    dpd_contract424(&E, &tIA, &ZMbEj, 3, 0, 0, -1, 0);
    dpd_buf4_close(&ZMbEj);
    dpd_buf4_close(&E);

    /* T1^2 -> Zmbej */

    /*** ABAB ***/
    dpd_buf4_init(&Y, CC_TMP0, 0, 0, 11, 0, 11, 0, "Y (Mn,Ej");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 10, 11, 10, 0, "CC2 ZMbEj (Ej,Mb)");
    dpd_contract424(&Y, &tIA, &Z, 1, 0, 0, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);

    /* Copy to Hbar */
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    dpd_buf4_sort_axpy(&Z, CC_TMP0, rspq, 11, 10, "CC2 ZMbEj (Ej,Mb)", 1);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ (Mb,Je)");
    dpd_buf4_sort(&Z, CC_TMP0, srqp, 11, 11, "CC2 ZMbeJ (eJ,bM)");
    dpd_buf4_close(&Z);
  }

  else if(params.ref == 1) { /** ROHF **/

    /* W(mb,je) <-- <mb||ej> */

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 ZMBEJ", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 Zmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 ZmBEj", -1);
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 ZMbeJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "CC2 ZMbEj");
    dpd_buf4_copy(&D, CC_TMP0, "CC2 ZmBeJ");
    dpd_buf4_close(&D);

    /* F -> Zmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&ZMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMBEJ");
    dpd_contract424(&F, &tIA, &ZMBEJ, 3, 1, 0, 1, 1);
    dpd_buf4_close(&ZMBEJ);
    dpd_buf4_init(&Zmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Zmbej");
    dpd_contract424(&F, &tia, &Zmbej, 3, 1, 0, 1, 1);
    dpd_buf4_close(&Zmbej);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&ZMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    dpd_contract424(&F, &tia, &ZMbEj, 3, 1, 0, 1, 1);
    dpd_buf4_close(&ZMbEj);
    dpd_buf4_init(&ZmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZmBeJ");
    dpd_contract424(&F, &tIA, &ZmBeJ, 3, 1, 0, 1, 1);
    dpd_buf4_close(&ZmBeJ);

    dpd_buf4_init(&ZMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
    dpd_contract244(&tIA, &F, &ZMbeJ, 1, 2, 1, -1, 1);
    dpd_buf4_close(&ZMbeJ);
    dpd_buf4_init(&ZmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
    dpd_contract244(&tia, &F, &ZmBEj, 1, 2, 1, -1, 1);
    dpd_buf4_close(&ZmBEj);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /* E -> Zmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    dpd_buf4_init(&ZMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMBEJ");
    dpd_contract424(&E, &tIA, &ZMBEJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&ZMBEJ);
    dpd_buf4_init(&Zmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Zmbej");
    dpd_contract424(&E, &tia, &Zmbej, 1, 0, 1, 1, 1);
    dpd_buf4_close(&Zmbej);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_init(&ZMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    dpd_contract424(&E, &tia, &ZMbEj, 3, 0, 1, -1, 1);
    dpd_buf4_close(&ZMbEj);
    dpd_buf4_init(&ZmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZmBeJ");
    dpd_contract424(&E, &tIA, &ZmBeJ, 3, 0, 1, -1, 1);
    dpd_buf4_close(&ZmBeJ);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_init(&ZMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
    dpd_contract424(&E, &tia, &ZMbeJ, 1, 0, 1, 1, 1);
    dpd_buf4_close(&ZMbeJ);
    dpd_buf4_init(&ZmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
    dpd_contract424(&E, &tIA, &ZmBEj, 1, 0, 1, 1, 1);
    dpd_buf4_close(&ZmBEj);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA);  
    dpd_file2_close(&tia);

    /* Convert to (ME,JB) for remaining terms */

    dpd_buf4_init(&ZMBEJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMBEJ");
    dpd_buf4_sort(&ZMBEJ, CC_TMP0, prsq, 10, 10, "CC2 ZMBEJ");
    dpd_buf4_close(&ZMBEJ);

    dpd_buf4_init(&Zmbej, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Zmbej");
    dpd_buf4_sort(&Zmbej, CC_TMP0, prsq, 10, 10, "CC2 Zmbej");
    dpd_buf4_close(&Zmbej);

    dpd_buf4_init(&ZMbEj, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    dpd_buf4_sort(&ZMbEj, CC_TMP0, prsq, 10, 10, "CC2 ZMbEj");
    dpd_buf4_close(&ZMbEj);

    dpd_buf4_init(&ZmBeJ, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZmBeJ");
    dpd_buf4_sort(&ZmBeJ, CC_TMP0, prsq, 10, 10, "CC2 ZmBeJ");
    dpd_buf4_close(&ZmBeJ);

    dpd_buf4_init(&ZMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
    dpd_buf4_sort(&ZMbeJ, CC_TMP0, psrq, 10, 10, "CC2 ZMbeJ");
    dpd_buf4_close(&ZMbeJ);

    dpd_buf4_init(&ZmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
    dpd_buf4_sort(&ZmBEj, CC_TMP0, psrq, 10, 10, "CC2 ZmBEj");
    dpd_buf4_close(&ZmBEj);

    /* X -> Zmbej */
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /*** AAAA ***/

    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMBEJ");
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
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMBEJ");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BBBB ***/

    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 Zmbej");
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
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 Zmbej");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABAB ***/
  
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbEj");
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
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbEj");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BABA ***/

    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBeJ");
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
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBeJ");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** ABBA ***/
  
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
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
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
    dpd_contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    /*** BAAB ***/

    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
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
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
    dpd_contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);


    /* D(me,nf) * T2(jb,nf) --> W(me,jb) */
    dpd_buf4_init(&Zmbej, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 Zmbej");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tiajb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &tiajb, &Zmbej, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiajb);

    /* D(me,NF) * T2(jb,NF) --> W(me,jb) */
    dpd_buf4_init(&tiaJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_contract444(&D, &tiaJB, &Zmbej, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiaJB);

    dpd_buf4_close(&Zmbej);


    /* D(ME,NF) * T2(JB,NF) --> W(ME,JB) */
    dpd_buf4_init(&ZMBEJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMBEJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tIAJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &tIAJB, &ZMBEJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAJB);

    /* D(ME,nf) * T2(JB,nf) --> W(ME,JB) */
    dpd_buf4_init(&tIAjb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_contract444(&D, &tIAjb, &ZMBEJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAjb);

    dpd_buf4_close(&ZMBEJ);


    /* D(me,nf) * T2(JB,nf) --> W(me,JB) */
    dpd_buf4_init(&ZmBeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBeJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tIAjb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &tIAjb, &ZmBeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAjb);

    /* D(me,NF) * T2(JB,NF) --> W(me,JB) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&tIAJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &tIAJB, &ZmBeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIAJB);

    dpd_buf4_close(&ZmBeJ);


    /* D(ME,NF) * T2(jb,NF) --> W(ME,jb) */
    dpd_buf4_init(&ZMbEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbEj");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&tiaJB, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&D, &tiaJB, &ZMbEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiaJB);

    /* D(ME,nf) * T2(jb,nf) --> W(ME,jb) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&tiajb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &tiajb, &ZMbEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tiajb);

    dpd_buf4_close(&ZMbEj);


    /* D(mE,Nf) * T2(jB,Nf) --> W(mE,jB) */
    dpd_buf4_init(&ZmBEj, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&tjAIb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_contract444(&D, &tjAIb, &ZmBEj, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tjAIb);

    dpd_buf4_close(&ZmBEj);


    /* D(Me,nF) * T2(Jb,nF) --> W(Me,Jb) */
    dpd_buf4_init(&ZMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&tIbjA, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_contract444(&D, &tIbjA, &ZMbeJ, 0, 0, 0.5, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIbjA);

    dpd_buf4_close(&ZMbeJ);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* C & D -> Zmbej */

    dpd_buf4_init(&C, CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 ZMBEJ", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 Zmbej", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_buf4_scmcopy(&D, CC_TMP0, "CC2 ZMbEj", -1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_buf4_scmcopy(&D, CC_TMP0, "CC2 ZmBeJ", -1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&C, CC_CINTS, 0, 27, 26, 27, 26, 0, "C <Ai|Bj> (iA,Bj)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 ZmBEj", -1);
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 24, 25, 24, 25, 0, "C <Ia|Jb> (Ia,bJ)");
    dpd_buf4_scmcopy(&C, CC_TMP0, "CC2 ZMbeJ", -1);
    dpd_buf4_close(&C);

    /* F -> Zmbej */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&Z, CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 ZMBEJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Zmbej");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_contract424(&F, &tia, &Z, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 ZMbEj");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract424(&F, &tia, &Z, 3, 1, 0, -1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 ZmBeJ");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, -1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Y, CC_TMP0, 0, 24, 24, 24, 24, 0, "CC2 YMbeJ (Mb,Je)");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract244(&tIA, &F, &Y, 1, 2, 1, -1, 0);
    dpd_buf4_close(&F);
    dpd_buf4_sort_axpy(&Y, CC_TMP0, pqsr, 24, 25, "CC2 ZMbeJ", 1);
    dpd_buf4_close(&Y);

    dpd_buf4_init(&Y, CC_TMP0, 0, 27, 27, 27, 27, 0, "CC2 YmBEj (mB,jE)");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract244(&tia, &F, &Y, 1, 2, 1, -1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_sort_axpy(&Y, CC_TMP0, pqsr, 27, 26, "CC2 ZmBEj", 1);
    dpd_buf4_close(&Y);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /* D & E -> Zmbej */

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

    dpd_buf4_init(&Z, CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 ZMBEJ");
    dpd_contract424(&Y, &tIA, &Z, 1, 0, 1, -0.5, 1);
    dpd_buf4_close(&Y);
    dpd_buf4_close(&Z);

    /*** BBBB ***/

    dpd_buf4_init(&E, CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    dpd_buf4_scmcopy(&E, CC_TMP0, "Ymnej", -1);
    dpd_buf4_close(&E);

    dpd_buf4_init(&Y, CC_TMP0, 0, 10, 31, 10, 31, 0, "Ymnej");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D, &tia, &Y, 3, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z, CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Zmbej");
    dpd_contract424(&Y, &tia, &Z, 1, 0, 1, -0.5, 1);
    dpd_buf4_close(&Y);
    dpd_buf4_close(&Z);

    /*** ABAB ***/

    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_buf4_copy(&E, CC_TMP0, "YMnEj");
    dpd_buf4_close(&E);

    dpd_buf4_init(&Y, CC_TMP0, 0, 22, 26, 22, 26, 0, "YMnEj");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &tia, &Y, 3, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z, CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 ZMbEj");
    dpd_contract424(&Y, &tia, &Z, 1, 0, 1, 0.5, 1);
    dpd_buf4_close(&Y);
    dpd_buf4_close(&Z);

    /*** BABA ***/

    dpd_buf4_init(&E, CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    dpd_buf4_copy(&E, CC_TMP0, "YmNeJ");
    dpd_buf4_close(&E);

    dpd_buf4_init(&Y, CC_TMP0, 0, 23, 25, 23, 25, 0, "YmNeJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract424(&D, &tIA, &Y, 3, 1, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z, CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 ZmBeJ");
    dpd_contract424(&Y, &tIA, &Z, 1, 0, 1, 0.5, 1);
    dpd_buf4_close(&Y);
    dpd_buf4_close(&Z);

    /*** ABBA ***/
  
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_copy(&E, CC_TMP0, "YMneJ (Mn,Je)");
    dpd_buf4_close(&E);

    dpd_buf4_init(&Y, CC_TMP0, 0, 22, 24, 22, 24, 0, "YMneJ (Mn,Je)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract244(&tIA, &D, &Y, 1, 2, 1, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Y1, CC_TMP0, 0, 24, 24, 24, 24, 0, "CC2 YMbeJ (Mb,Je)");
    dpd_contract424(&Y, &tia, &Y1, 1, 0, 1, 1, 0);
    dpd_buf4_close(&Y);
    dpd_buf4_sort_axpy(&Y1, CC_TMP0, pqsr, 24, 25, "CC2 ZMbeJ", 0.5);
    dpd_buf4_close(&Y1);

    /*** BAAB ***/
  
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_buf4_copy(&E, CC_TMP0, "YmNEj (mN,jE)");
    dpd_buf4_close(&E);

    dpd_buf4_init(&Y, CC_TMP0, 0, 23, 27, 23, 27, 0, "YmNEj (mN,jE)");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract244(&tia, &D, &Y, 1, 2, 1, 1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Y1, CC_TMP0, 0, 27, 27, 27, 27, 0, "CC2 YmBEj (mB,jE)");
    dpd_contract424(&Y, &tIA, &Y1, 1, 0, 1, 1, 0);
    dpd_buf4_close(&Y);
    dpd_buf4_sort_axpy(&Y1, CC_TMP0, pqsr, 27, 26, "CC2 ZmBEj", 0.5);
    dpd_buf4_close(&Y1);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

  } /** UHF **/

  return;
}

}} // namespace psi::cchbar
