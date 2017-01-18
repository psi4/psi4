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
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
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

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "CC2 ZMbEj");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 ZMbeJ (Mb,Je)", 1);
    global_dpd_->buf4_close(&C);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /* F -> Zmbej */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&ZMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    global_dpd_->contract424(&F, &tIA, &ZMbEj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&ZMbEj);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    global_dpd_->contract424(&F, &tIA, &Z, 3, 1, 0, -1, 0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 10, 10, "Z(Mb,Je)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Je)");
    global_dpd_->buf4_init(&ZMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ (Mb,Je)");
    global_dpd_->buf4_axpy(&Z, &ZMbeJ, -1.0);
    global_dpd_->buf4_close(&ZMbeJ);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&F);

    /* E -> Zmbej */

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&ZMbEj, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "CC2 ZMbEj (Ej,Mb)");
    global_dpd_->contract424(&E, &tIA, &ZMbEj, 3, 0, 0, -1, 0);
    global_dpd_->buf4_close(&ZMbEj);
    global_dpd_->buf4_close(&E);

    /* T1^2 -> Zmbej */

    /*** ABAB ***/
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Y (Mn,Ej");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "CC2 ZMbEj (Ej,Mb)");
    global_dpd_->contract424(&Y, &tIA, &Z, 1, 0, 0, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);

    /* Copy to Hbar */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, rspq, 11, 10, "CC2 ZMbEj (Ej,Mb)", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ (Mb,Je)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, srqp, 11, 11, "CC2 ZMbeJ (eJ,bM)");
    global_dpd_->buf4_close(&Z);
  }

  else if(params.ref == 1) { /** ROHF **/

    /* W(mb,je) <-- <mb||ej> */

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 ZMBEJ", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 Zmbej", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 ZmBEj", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 ZMbeJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "CC2 ZMbEj");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "CC2 ZmBeJ");
    global_dpd_->buf4_close(&D);

    /* F -> Zmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&ZMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMBEJ");
    global_dpd_->contract424(&F, &tIA, &ZMBEJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&ZMBEJ);
    global_dpd_->buf4_init(&Zmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Zmbej");
    global_dpd_->contract424(&F, &tia, &Zmbej, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Zmbej);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&ZMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    global_dpd_->contract424(&F, &tia, &ZMbEj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&ZMbEj);
    global_dpd_->buf4_init(&ZmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZmBeJ");
    global_dpd_->contract424(&F, &tIA, &ZmBeJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&ZmBeJ);

    global_dpd_->buf4_init(&ZMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
    global_dpd_->contract244(&tIA, &F, &ZMbeJ, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&ZMbeJ);
    global_dpd_->buf4_init(&ZmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
    global_dpd_->contract244(&tia, &F, &ZmBEj, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&ZmBEj);
    global_dpd_->buf4_close(&F);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* E -> Zmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_init(&ZMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMBEJ");
    global_dpd_->contract424(&E, &tIA, &ZMBEJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&ZMBEJ);
    global_dpd_->buf4_init(&Zmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Zmbej");
    global_dpd_->contract424(&E, &tia, &Zmbej, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Zmbej);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&ZMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    global_dpd_->contract424(&E, &tia, &ZMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&ZMbEj);
    global_dpd_->buf4_init(&ZmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZmBeJ");
    global_dpd_->contract424(&E, &tIA, &ZmBeJ, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&ZmBeJ);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&ZMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
    global_dpd_->contract424(&E, &tia, &ZMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&ZMbeJ);
    global_dpd_->buf4_init(&ZmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
    global_dpd_->contract424(&E, &tIA, &ZmBEj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&ZmBEj);
    global_dpd_->buf4_close(&E);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* Convert to (ME,JB) for remaining terms */

    global_dpd_->buf4_init(&ZMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMBEJ");
    global_dpd_->buf4_sort(&ZMBEJ, PSIF_CC_TMP0, prsq, 10, 10, "CC2 ZMBEJ");
    global_dpd_->buf4_close(&ZMBEJ);

    global_dpd_->buf4_init(&Zmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Zmbej");
    global_dpd_->buf4_sort(&Zmbej, PSIF_CC_TMP0, prsq, 10, 10, "CC2 Zmbej");
    global_dpd_->buf4_close(&Zmbej);

    global_dpd_->buf4_init(&ZMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZMbEj");
    global_dpd_->buf4_sort(&ZMbEj, PSIF_CC_TMP0, prsq, 10, 10, "CC2 ZMbEj");
    global_dpd_->buf4_close(&ZMbEj);

    global_dpd_->buf4_init(&ZmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 ZmBeJ");
    global_dpd_->buf4_sort(&ZmBeJ, PSIF_CC_TMP0, prsq, 10, 10, "CC2 ZmBeJ");
    global_dpd_->buf4_close(&ZmBeJ);

    global_dpd_->buf4_init(&ZMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
    global_dpd_->buf4_sort(&ZMbeJ, PSIF_CC_TMP0, psrq, 10, 10, "CC2 ZMbeJ");
    global_dpd_->buf4_close(&ZMbeJ);

    global_dpd_->buf4_init(&ZmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
    global_dpd_->buf4_sort(&ZmBEj, PSIF_CC_TMP0, psrq, 10, 10, "CC2 ZmBEj");
    global_dpd_->buf4_close(&ZmBEj);

    /* X -> Zmbej */
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /*** AAAA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMBEJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMBEJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BBBB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 Zmbej");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 Zmbej");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABAB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbEj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbEj");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BABA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBeJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBeJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABBA ***/

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BAAB ***/

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
    global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract444(&D, &t2, &W, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);


    /* D(me,nf) * T2(jb,nf) --> W(me,jb) */
    global_dpd_->buf4_init(&Zmbej, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 Zmbej");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&tiajb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &tiajb, &Zmbej, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tiajb);

    /* D(me,NF) * T2(jb,NF) --> W(me,jb) */
    global_dpd_->buf4_init(&tiaJB, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->contract444(&D, &tiaJB, &Zmbej, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tiaJB);

    global_dpd_->buf4_close(&Zmbej);


    /* D(ME,NF) * T2(JB,NF) --> W(ME,JB) */
    global_dpd_->buf4_init(&ZMBEJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMBEJ");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&tIAJB, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &tIAJB, &ZMBEJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIAJB);

    /* D(ME,nf) * T2(JB,nf) --> W(ME,JB) */
    global_dpd_->buf4_init(&tIAjb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->contract444(&D, &tIAjb, &ZMBEJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIAjb);

    global_dpd_->buf4_close(&ZMBEJ);


    /* D(me,nf) * T2(JB,nf) --> W(me,JB) */
    global_dpd_->buf4_init(&ZmBeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBeJ");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&tIAjb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &tIAjb, &ZmBeJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIAjb);

    /* D(me,NF) * T2(JB,NF) --> W(me,JB) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&tIAJB, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &tIAJB, &ZmBeJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIAJB);

    global_dpd_->buf4_close(&ZmBeJ);


    /* D(ME,NF) * T2(jb,NF) --> W(ME,jb) */
    global_dpd_->buf4_init(&ZMbEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbEj");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&tiaJB, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &tiaJB, &ZMbEj, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tiaJB);

    /* D(ME,nf) * T2(jb,nf) --> W(ME,jb) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&tiajb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &tiajb, &ZMbEj, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tiajb);

    global_dpd_->buf4_close(&ZMbEj);


    /* D(mE,Nf) * T2(jB,Nf) --> W(mE,jB) */
    global_dpd_->buf4_init(&ZmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZmBEj");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_init(&tjAIb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->contract444(&D, &tjAIb, &ZmBEj, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tjAIb);

    global_dpd_->buf4_close(&ZmBEj);


    /* D(Me,nF) * T2(Jb,nF) --> W(Me,Jb) */
    global_dpd_->buf4_init(&ZMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 ZMbeJ");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_init(&tIbjA, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->contract444(&D, &tIbjA, &ZMbeJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIbjA);

    global_dpd_->buf4_close(&ZMbeJ);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* C & D -> Zmbej */

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 ZMBEJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 Zmbej", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "CC2 ZMbEj", -1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "CC2 ZmBeJ", -1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 26, 27, 26, 0, "C <Ai|Bj> (iA,Bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 ZmBEj", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 25, 24, 25, 0, "C <Ia|Jb> (Ia,bJ)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 ZMbeJ", -1);
    global_dpd_->buf4_close(&C);

    /* F -> Zmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 ZMBEJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract424(&F, &tIA, &Z, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Zmbej");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &tia, &Z, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 ZMbEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract424(&F, &tia, &Z, 3, 1, 0, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 ZmBeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract424(&F, &tIA, &Z, 3, 1, 0, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "CC2 YMbeJ (Mb,Je)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract244(&tIA, &F, &Y, 1, 2, 1, -1, 0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_sort_axpy(&Y, PSIF_CC_TMP0, pqsr, 24, 25, "CC2 ZMbeJ", 1);
    global_dpd_->buf4_close(&Y);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "CC2 YmBEj (mB,jE)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract244(&tia, &F, &Y, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_sort_axpy(&Y, PSIF_CC_TMP0, pqsr, 27, 26, "CC2 ZmBEj", 1);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* D & E -> Zmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /*** AAAA ***/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    global_dpd_->buf4_scmcopy(&E, PSIF_CC_TMP0, "YMNEJ", -1);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 0, 21, 0, 21, 0, "YMNEJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 ZMBEJ");
    global_dpd_->contract424(&Y, &tIA, &Z, 1, 0, 1, -0.5, 1);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_close(&Z);

    /*** BBBB ***/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_scmcopy(&E, PSIF_CC_TMP0, "Ymnej", -1);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 31, 10, 31, 0, "Ymnej");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D, &tia, &Y, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Zmbej");
    global_dpd_->contract424(&Y, &tia, &Z, 1, 0, 1, -0.5, 1);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_close(&Z);

    /*** ABAB ***/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->buf4_copy(&E, PSIF_CC_TMP0, "YMnEj");
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 22, 26, 22, 26, 0, "YMnEj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract424(&D, &tia, &Y, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 ZMbEj");
    global_dpd_->contract424(&Y, &tia, &Z, 1, 0, 1, 0.5, 1);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_close(&Z);

    /*** BABA ***/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    global_dpd_->buf4_copy(&E, PSIF_CC_TMP0, "YmNeJ");
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 23, 25, 23, 25, 0, "YmNeJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 ZmBeJ");
    global_dpd_->contract424(&Y, &tIA, &Z, 1, 0, 1, 0.5, 1);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_close(&Z);

    /*** ABBA ***/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_copy(&E, PSIF_CC_TMP0, "YMneJ (Mn,Je)");
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 22, 24, 22, 24, 0, "YMneJ (Mn,Je)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Y1, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "CC2 YMbeJ (Mb,Je)");
    global_dpd_->contract424(&Y, &tia, &Y1, 1, 0, 1, 1, 0);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_sort_axpy(&Y1, PSIF_CC_TMP0, pqsr, 24, 25, "CC2 ZMbeJ", 0.5);
    global_dpd_->buf4_close(&Y1);

    /*** BAAB ***/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->buf4_copy(&E, PSIF_CC_TMP0, "YmNEj (mN,jE)");
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 23, 27, 23, 27, 0, "YmNEj (mN,jE)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Y1, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "CC2 YmBEj (mB,jE)");
    global_dpd_->contract424(&Y, &tIA, &Y1, 1, 0, 1, 1, 0);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_sort_axpy(&Y1, PSIF_CC_TMP0, pqsr, 27, 26, "CC2 ZmBEj", 0.5);
    global_dpd_->buf4_close(&Y1);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  } /** UHF **/

  return;
}

}} // namespace psi::cchbar
