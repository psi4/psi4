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

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 WMbeJ (Mb,Je)", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "CC2 WMbEj");
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /* F -> Wmbej */
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    global_dpd_->contract424(&F, &tIA, &WMbEj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    global_dpd_->contract424(&F, &tIA, &Z, 3, 1, 0, -1, 0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 10, 10, "Z(Mb,Je)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Je)");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Mb,Je)");
    global_dpd_->buf4_axpy(&Z, &WMbeJ, 1.0);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&F);

    /* E -> Wmbej */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    global_dpd_->contract424(&E, &tIA, &WMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Mb,Je)");
    global_dpd_->contract424(&E, &tIA, &WMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_close(&E);

    global_dpd_->file2_close(&tIA);


    /* Sort to (ME,JB) */
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC2_HET1, prsq, 10, 10, "CC2 WMbEj (ME,jb)");
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Mb,Je)");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC2_HET1, psrq, 10, 10, "CC2 WMbeJ (Me,Jb)");
    global_dpd_->buf4_close(&WMbeJ);


    /* T1^2 -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /*** ABAB ***/
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj (ME,jb)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABBA ***/
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Me,Jb)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);

  }

  else if(params.ref == 1) { /** ROHF **/

    /* W(mb,je) <-- <mb||ej> */

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 WMBEJ", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 Wmbej", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 WmBEj", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 WMbeJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "CC2 WMbEj");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "CC2 WmBeJ");
    global_dpd_->buf4_close(&D);

    /* F -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMBEJ");
    global_dpd_->contract424(&F, &tIA, &WMBEJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Wmbej");
    global_dpd_->contract424(&F, &tia, &Wmbej, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    global_dpd_->contract424(&F, &tia, &WMbEj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WmBeJ");
    global_dpd_->contract424(&F, &tIA, &WmBeJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WmBeJ);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    global_dpd_->contract244(&tIA, &F, &WMbeJ, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
    global_dpd_->contract244(&tia, &F, &WmBEj, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&F);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* E -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMBEJ");
    global_dpd_->contract424(&E, &tIA, &WMBEJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Wmbej");
    global_dpd_->contract424(&E, &tia, &Wmbej, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    global_dpd_->contract424(&E, &tia, &WMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WmBeJ");
    global_dpd_->contract424(&E, &tIA, &WmBeJ, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WmBeJ);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    global_dpd_->contract424(&E, &tia, &WMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
    global_dpd_->contract424(&E, &tIA, &WmBEj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&E);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* Convert to (ME,JB) for remaining terms */

    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMBEJ");
    global_dpd_->buf4_sort(&WMBEJ, PSIF_CC2_HET1, prsq, 10, 10, "CC2 WMBEJ");
    global_dpd_->buf4_close(&WMBEJ);

    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 Wmbej");
    global_dpd_->buf4_sort(&Wmbej, PSIF_CC2_HET1, prsq, 10, 10, "CC2 Wmbej");
    global_dpd_->buf4_close(&Wmbej);

    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC2_HET1, prsq, 10, 10, "CC2 WMbEj");
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "CC2 WmBeJ");
    global_dpd_->buf4_sort(&WmBeJ, PSIF_CC2_HET1, prsq, 10, 10, "CC2 WmBeJ");
    global_dpd_->buf4_close(&WmBeJ);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC2_HET1, psrq, 10, 10, "CC2 WMbeJ");
    global_dpd_->buf4_close(&WMbeJ);

    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
    global_dpd_->buf4_sort(&WmBEj, PSIF_CC2_HET1, psrq, 10, 10, "CC2 WmBEj");
    global_dpd_->buf4_close(&WmBEj);

    /* X -> Wmbej */
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /*** AAAA ***/

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMBEJ");
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
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMBEJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BBBB ***/

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 Wmbej");
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
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 Wmbej");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABAB ***/

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj");
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
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BABA ***/

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBeJ");
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
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBeJ");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABBA ***/

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
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
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BAAB ***/

    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
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
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);


    /* D(me,nf) * T2(jb,nf) --> W(me,jb) */
    global_dpd_->buf4_init(&Wmbej, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 Wmbej");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&tiajb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &tiajb, &Wmbej, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tiajb);

    /* D(me,NF) * T2(jb,NF) --> W(me,jb) */
    global_dpd_->buf4_init(&tiaJB, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->contract444(&D, &tiaJB, &Wmbej, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tiaJB);

    global_dpd_->buf4_close(&Wmbej);


    /* D(ME,NF) * T2(JB,NF) --> W(ME,JB) */
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMBEJ");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&tIAJB, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &tIAJB, &WMBEJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIAJB);

    /* D(ME,nf) * T2(JB,nf) --> W(ME,JB) */
    global_dpd_->buf4_init(&tIAjb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->contract444(&D, &tIAjb, &WMBEJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIAjb);

    global_dpd_->buf4_close(&WMBEJ);


    /* D(me,nf) * T2(JB,nf) --> W(me,JB) */
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBeJ");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&tIAjb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->contract444(&D, &tIAjb, &WmBeJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIAjb);

    /* D(me,NF) * T2(JB,NF) --> W(me,JB) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&tIAJB, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->contract444(&D, &tIAJB, &WmBeJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIAJB);

    global_dpd_->buf4_close(&WmBeJ);


    /* D(ME,NF) * T2(jb,NF) --> W(ME,jb) */
    global_dpd_->buf4_init(&WMbEj, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_init(&tiaJB, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->contract444(&D, &tiaJB, &WMbEj, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tiaJB);

    /* D(ME,nf) * T2(jb,nf) --> W(ME,jb) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_init(&tiajb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->contract444(&D, &tiajb, &WMbEj, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tiajb);

    global_dpd_->buf4_close(&WMbEj);


    /* D(mE,Nf) * T2(jB,Nf) --> W(mE,jB) */
    global_dpd_->buf4_init(&WmBEj, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WmBEj");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_init(&tjAIb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->contract444(&D, &tjAIb, &WmBEj, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tjAIb);

    global_dpd_->buf4_close(&WmBEj);


    /* D(Me,nF) * T2(Jb,nF) --> W(Me,Jb) */
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_init(&tIbjA, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->contract444(&D, &tIbjA, &WMbeJ, 0, 0, 0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&tIbjA);

    global_dpd_->buf4_close(&WMbeJ);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* W(mb,je) <-- <mb||ej> */

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 WMBEJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "CC2 Wmbej", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "CC2 WMbEj", -1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "CC2 WmBeJ", -1);
    global_dpd_->buf4_close(&D);

    /* F -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 WMBEJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Wmbej");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 WMbEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract424(&F, &tia, &W, 3, 1, 0, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 WmBeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract424(&F, &tIA, &W, 3, 1, 0, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* D & E -> Wmbej */

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

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 WMBEJ");
    global_dpd_->contract424(&Y, &tIA, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_close(&W);

    /*** BBBB ***/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_scmcopy(&E, PSIF_CC_TMP0, "Ymnej", -1);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 31, 10, 31, 0, "Ymnej");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D, &tia, &Y, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Wmbej");
    global_dpd_->contract424(&Y, &tia, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_close(&W);

    /*** ABAB ***/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->buf4_copy(&E, PSIF_CC_TMP0, "YMnEj");
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 22, 26, 22, 26, 0, "YMnEj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract424(&D, &tia, &Y, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 WMbEj");
    global_dpd_->contract424(&Y, &tia, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_close(&W);

    /*** BABA ***/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    global_dpd_->buf4_copy(&E, PSIF_CC_TMP0, "YmNeJ");
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 23, 25, 23, 25, 0, "YmNeJ");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 WmBeJ");
    global_dpd_->contract424(&Y, &tIA, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Y);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /** Sort **/
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "CC2 WMBEJ");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, prqs, 20, 21, "CC2 WMBEJ (ME,BJ)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "CC2 Wmbej");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, prqs, 30, 31, "CC2 Wmbej (me,bj)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "CC2 WMbEj");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, prqs, 20, 31, "CC2 WMbEj (ME,bj)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "CC2 WmBeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC2_HET1, prqs, 30, 21, "CC2 WmBeJ (me,BJ)");
    global_dpd_->buf4_close(&W);

  } /** UHF **/

  return;
}

}} // namespace psi::cchbar
