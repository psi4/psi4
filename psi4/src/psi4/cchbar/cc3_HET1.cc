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

/* compute remaining matrix elements for [ H, e^T1 ] which are not needed by
 * ccenergy for CC3 but are needed for CC3 EOM and CC3 Lambda and CC3 response.
 *
 * Wmbej = <mb||ej> + t_j^f <mb||ef> - t_n^b <mn||ej> - t_j^f t_n^b <mn||ef>
 *
*/

void HET1_Wmbej(void);
void purge_HET1_Wmbej(void);
void HET1_Wabef(void);

void cc3_HET1(void) {

  HET1_Wmbej();

  if (params.ref == 1) purge_HET1_Wmbej();

  HET1_Wabef();

  return;
}


void HET1_Wmbej(void)
{
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ;
  dpdbuf4 D, C, F, E, X, Y, t2, W, Z;
  dpdfile2 tIA, tia;
  int i;

  if(params.ref == 0) { /** RHF **/

    /* <mb||ej> -> Wmbej */

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMbeJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "WMbEj");
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /* F -> Wmbej */

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");

    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&F, &tIA, &WMbEj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(bM,eJ)");
    global_dpd_->contract424(&F, &tIA, &Z, 3, 1, 0, -1, 0);
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 10, 10, "Z(Mb,Je)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Je)");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_axpy(&Z, &WMbeJ, 1.0);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_close(&F);

    global_dpd_->file2_close(&tIA);

    /* E -> Wmbej */
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&E, &tIA, &WMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&E, &tIA, &WMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_close(&E);

    global_dpd_->file2_close(&tIA);


    /* Sort to (ME,JB) */

    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC3_HET1, prsq, 10, 10, "CC3 WMbEj (ME,jb)");
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC3_HET1, psrq, 10, 10, "CC3 WMbeJ (Me,Jb)");
    global_dpd_->buf4_close(&WMbeJ);


    /* T1^2 -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);

    /* also store WMbEj (Mb,Ej)  WMbeJ (bM,eJ)*/
    global_dpd_->buf4_init(&WMbEj, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC3_HET1, psqr, 10, 11, "CC3 WMbEj (Mb,Ej)");
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WMbEj, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC3_HET1, spqr, 11, 11, "CC3 WMbeJ (bM,eJ)");
    global_dpd_->buf4_close(&WMbEj);

  }
  else if(params.ref == 1) { /** ROHF **/

    /* W(mb,je) <-- <mb||ej> */

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMBEJ", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "Wmbej", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WmBEj", -1);
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMbeJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "WMbEj");
    global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "WmBeJ");
    global_dpd_->buf4_close(&D);

    /* F -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->contract424(&F, &tIA, &WMBEJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->contract424(&F, &tia, &Wmbej, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&F);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&F, &tia, &WMbEj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->contract424(&F, &tIA, &WmBeJ, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&WmBeJ);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract244(&tIA, &F, &WMbeJ, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract244(&tia, &F, &WmBEj, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&F);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* E -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->contract424(&E, &tIA, &WMBEJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMBEJ);
    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->contract424(&E, &tia, &Wmbej, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Wmbej);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->contract424(&E, &tia, &WMbEj, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WMbEj);
    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->contract424(&E, &tIA, &WmBeJ, 3, 0, 1, -1, 1);
    global_dpd_->buf4_close(&WmBeJ);
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract424(&E, &tia, &WMbeJ, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMbeJ);
    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract424(&E, &tIA, &WmBEj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WmBEj);
    global_dpd_->buf4_close(&E);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* Convert to (ME,JB) for remaining terms */

    global_dpd_->buf4_init(&WMBEJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMBEJ");
    global_dpd_->buf4_sort(&WMBEJ, PSIF_CC3_HET1, prsq, 10, 10, "CC3 WMBEJ (ME,JB)");
    global_dpd_->buf4_close(&WMBEJ);

    global_dpd_->buf4_init(&Wmbej, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Wmbej");
    global_dpd_->buf4_sort(&Wmbej, PSIF_CC3_HET1, prsq, 10, 10, "CC3 Wmbej (me,jb)");
    global_dpd_->buf4_close(&Wmbej);

    global_dpd_->buf4_init(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj");
    global_dpd_->buf4_sort(&WMbEj, PSIF_CC3_HET1, prsq, 10, 10, "CC3 WMbEj (ME,jb)");
    global_dpd_->buf4_close(&WMbEj);

    global_dpd_->buf4_init(&WmBeJ, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WmBeJ");
    global_dpd_->buf4_sort(&WmBeJ, PSIF_CC3_HET1, prsq, 10, 10, "CC3 WmBeJ (me,JB)");
    global_dpd_->buf4_close(&WmBeJ);

    global_dpd_->buf4_init(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&WMbeJ, PSIF_CC3_HET1, psrq, 10, 10, "CC3 WMbeJ (Me,Jb)");
    global_dpd_->buf4_close(&WMbeJ);

    global_dpd_->buf4_init(&WmBEj, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->buf4_sort(&WmBEj, PSIF_CC3_HET1, psrq, 10, 10, "CC3 WmBEj (mE,jB)");
    global_dpd_->buf4_close(&WmBEj);

    /* X -> Wmbej */
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /*** AAAA ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMBEJ (ME,JB)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 10, 11, "CC3 WMBEJ (MB,EJ)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BBBB ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 Wmbej (me,jb)");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 10, 11, "CC3 Wmbej (mb,ej)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABAB ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 10, 11, "CC3 WMbEj (Mb,Ej)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BABA ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WmBeJ (me,JB)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 10, 11, "CC3 WmBeJ (mB,eJ)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABBA ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 10, 11, "CC3 WMbeJ (Mb,eJ)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BAAB ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WmBEj (mE,jB)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 10, 11, "CC3 WmBEj (mB,Ej)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  } /** ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* W(mb,je) <-- <mb||ej> */

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMBEJ", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "Wmbej", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "WMbEj", 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_TMP0, "WmBeJ", 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WmBEj", -1);
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    global_dpd_->buf4_scmcopy(&C, PSIF_CC_TMP0, "WMbeJ", -1);
    global_dpd_->buf4_close(&C);

    /* F -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    global_dpd_->contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    global_dpd_->contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract424(&F, &tia, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract424(&F, &tIA, &W, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    global_dpd_->contract244(&tIA, &F, &W, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    global_dpd_->contract244(&tia, &F, &W, 1, 2, 1, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* E -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    global_dpd_->contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->contract424(&E, &tia, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    global_dpd_->contract424(&E, &tIA, &W, 1, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->contract424(&E, &tia, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->contract424(&E, &tIA, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* Convert to (ME,JB) for remaining terms */

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "WMBEJ");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, prsq, 20, 20, "CC3 WMBEJ (ME,JB)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Wmbej");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, prsq, 30, 30, "CC3 Wmbej (me,jb)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "WMbEj");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, prsq, 20, 30, "CC3 WMbEj (ME,jb)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "WmBeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, prsq, 30, 20, "CC3 WmBeJ (me,JB)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psrq, 24, 24, "CC3 WMbeJ (Me,Jb)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psrq, 27, 27, "CC3 WmBEj (mE,jB)");
    global_dpd_->buf4_close(&W);

    /* X -> Wmbej */

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /*** AAAA ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Y (ME,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 21, 20, 21, 0, "D <IJ||AB> (IA,BJ)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 20, 20, 20, 0, "CC3 WMBEJ (ME,JB)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BBBB ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Y (me,jn)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 31, 30, 31, 0, "D <ij||ab> (ia,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 30, 30, 30, 0, "CC3 Wmbej (me,jb)");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABAB ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 20, 10, 20, 10, 0, "Y (ME,jn)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 31, 20, 31, 0, "D <Ij|Ab> (IA,bj)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 30, 20, 30, 0, "CC3 WMbEj (ME,jb)");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BABA ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 30, 0, 30, 0, 0, "Y (me,JN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 21, 30, 21, 0, "D <Ij|Ab> (ia,BJ)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 20, 30, 20, 0, "CC3 WmBeJ (me,JB)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** ABBA ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Y (Me,Jn)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->contract244(&tIA, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 24, 24, 24, 24, 0, "CC3 WMbeJ (Me,Jb)");
    global_dpd_->contract424(&Y, &tia, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    /*** BAAB ***/

    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "Y (mE,jN)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->contract244(&tia, &D, &Y, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 27, 27, 27, 27, 0, "CC3 WmBEj (mE,jB)");
    global_dpd_->contract424(&Y, &tIA, &W, 3, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* also store lists as Wmbej (mb,ej) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 20, 20, 20, 0, "CC3 WMBEJ (ME,JB)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 20, 21, "CC3 WMBEJ (MB,EJ)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 30, 30, 30, 0, "CC3 Wmbej (me,jb)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 30, 31, "CC3 Wmbej (mb,ej)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 20, 30, 20, 30, 0, "CC3 WMbEj (ME,jb)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 24, 26, "CC3 WMbEj (Mb,Ej)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 30, 20, 30, 20, 0, "CC3 WmBeJ (me,JB)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 27, 25, "CC3 WmBeJ (mB,eJ)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 24, 24, 24, 24, 0, "CC3 WMbeJ (Me,Jb)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 24, 25, "CC3 WMbeJ (Mb,eJ)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 27, 27, 27, 27, 0, "CC3 WmBEj (mE,jB)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psqr, 27, 26, "CC3 WmBEj (mB,Ej)");
    global_dpd_->buf4_close(&W);

    /* also make some Wmbej (bm,ej) */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 24, 25, 24, 25, 0, "CC3 WMbeJ (Mb,eJ)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qprs, 25, 25, "CC3 WMbeJ (bM,eJ)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 27, 26, 27, 26, 0, "CC3 WmBEj (mB,Ej)");
    global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, qprs, 26, 26, "CC3 WmBEj (Bm,Ej)");
    global_dpd_->buf4_close(&W);

  } /** UHF **/
  return;
}



/* Purge Wmbej matrix elements for ROHF */
void purge_HET1_Wmbej(void) {
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

  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 10, 10,"CC3 WMBEJ (ME,JB)");
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


  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 10, 10,"CC3 Wmbej (me,jb)");
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


  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 10, 10,"CC3 WMbEj (ME,jb)");
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


  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 10, 10,"CC3 WmBeJ (me,JB)");
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


  global_dpd_->file4_init(&W, PSIF_CC3_HET1, 0, 10, 10,"CC3 WmBEj (mE,jB)");
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
  return;
}

}} // namespace psi::cchbar
