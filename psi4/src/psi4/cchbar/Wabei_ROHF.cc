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

/** Wabei intermediates are stored here as (ei,ab) **/

void Wabei_ROHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  /** <EI||AB> **/
  global_dpd_->buf4_sort(&F, PSIF_CC_HBAR, qprs, 11, 7, "WEIAB");
  /** <ei||ab> **/
  global_dpd_->buf4_sort(&F, PSIF_CC_HBAR, qprs, 11, 7, "Weiab");
  global_dpd_->buf4_close(&F);

  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WEIAB");
  global_dpd_->buf4_scm(&W, -1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Weiab");
  global_dpd_->buf4_scm(&W, -1.0);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  /** <iE|bA> **/
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, qprs, 11, 5, "W(Ei,bA)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "W(Ei,bA)");
  global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WEiAb");
  global_dpd_->buf4_close(&W);
  /** <Ie|Ba> **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  global_dpd_->buf4_copy(&W, PSIF_CC_HBAR, "WeIaB");
  global_dpd_->buf4_close(&W);

  /** - F_ME t_MI^AB **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WEIAB");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /** - F_me t_mi^ab **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Weiab");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /** - F_ME t_Mi^Ab **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /** - F_me t_mI^aB **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /** <AB||EF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 7, 11, 7, 11, 0, "Z(AB,EI)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <ab|cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, rspq, 11, 7, "Z(EI,AB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 7, 11, 7, 0, "Z(EI,AB)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WEIAB");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /** <ab||ef> t_i^f **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 7, 11, 7, 11, 0, "Z(ab,ei)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <ab|cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, rspq, 11, 7, "Z(ei,ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 7, 11, 7, 0, "Z(ei,ab)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Weiab");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /** <Ab|Ef> t_i^f **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 11, 5, 11, 0, "Z(Ab,Ei)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, rspq, 11, 5, "Z(Ei,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(Ei,Ab)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /** - <Fe|Ba> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 11, 5, 11, 0, "Z(aB,eI)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, rspq, 11, 5, "Z(eI,aB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(eI,aB)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /** Prepare intermediates for second Wabef contribution to Wabei **/

  /** Z(MA,EI) <-- <MA||EF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(MA,EI)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** t_M^B Z(MA,EI) --> Z'(BA,EI) --> Z1(AB,EI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 5, 11, 5, 11, 0, "Z(BA,EI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(MA,EI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP2, qprs, 5, 11, "Z(AB,EI)");
  global_dpd_->buf4_close(&Z1);

  /** t_M^A Z(MB,EI) --> Z2(AB,EI) **/
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 5, 11, 5, 11, 0, "Z(AB,EI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(MA,EI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);

  /** Z1(AB,EI) + Z2(AB,EI) --> W(AB,EI) --> W(EI,AB) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z(AB,EI)");
  global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, rspq, 11, 5, "Z(EI,AB)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(EI,AB)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "WEIAB");
  global_dpd_->buf4_axpy(&Z2, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z2);


  /** Z(ma,ei) <-- <ma||ef> t_i^f **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(ma,ei)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** t_m^b Z(ma,ei) --> Z'(ba,ei) --> Z1(ab,ei) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 5, 11, 5, 11, 0, "Z(ba,ei)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(ma,ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP2, qprs, 5, 11, "Z(ab,ei)");
  global_dpd_->buf4_close(&Z1);

  /** t_m^a Z(mb,ei) --> Z2(ab,ei) **/
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 5, 11, 5, 11, 0, "Z(ab,ei)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(ma,ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);

  /** Z1(ab,ei) + Z2(ab,ei) --> W(ab,ei) --> W(ie,ab) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z(ab,ei)");
  global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, rspq, 11, 5, "Z(ei,ab)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(ei,ab)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "Weiab");
  global_dpd_->buf4_axpy(&Z2, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z2);

  /** t_i^f <mA|fE> --> Z(mA,Ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(mA,iE)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract244(&T1, &F, &Z, 1, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, pqsr, 10, 11, "Z(mA,Ei)");
  global_dpd_->buf4_close(&Z);

  /** <Mb|Ef> t_i^f --> Z(Mb,Ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(Mb,Ei)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** - T_M^A Z(Mb,Ei) --> Z(Ab,Ei) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z(Ab,Ei)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(Mb,Ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&Z1);

  /** - t_m^b Z(mA,Ei) --> Z1(Ab,Ei) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 11, 5, 11, 0, "Z1(bA,Ei)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(mA,Ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qprs, 5, 11, "Z(Ab,Ei)");
  global_dpd_->buf4_close(&Z1);

  /** Z1(Ab,Ei) + Z2(Ab,Ei) --> W(Ab,Ei) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 5, 11, 5, 11, 0, "Z(Ab,Ei)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z(Ab,Ei)");
  global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, rspq, 11, 5, "Z(Ei,Ab)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(Ei,Ab)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** t_I^F <Ma|Fe> --> Z(Ma,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Ma,Ie)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &F, &Z, 1, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, pqsr, 10, 11, "Z(Ma,eI)");
  global_dpd_->buf4_close(&Z);

  /** <mB|eF> t_I^F --> Z(mB,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(mB,eI)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** t_m^a Z(mB,eI) --> Z(aB,eI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z(aB,eI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(mB,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&Z1);

  /** t_M^B Z(Ma,eI) --> Z1(aB,eI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 11, 5, 11, 0, "Z1(Ba,eI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "Z(Ma,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qprs, 5, 11, "Z(aB,eI)");
  global_dpd_->buf4_close(&Z1);

  /** Z1(aB,eI) + Z2(aB,eI) --> W(aB,eI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 5, 11, 5, 11, 0, "Z(aB,eI)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z(aB,eI)");
  global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, rspq, 11, 5, "Z(eI,aB)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(eI,aB)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);



  /** Final term of Wabef contribution to Wabei **/

  /** t_I^F <MN||EF> --> Z(MN,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 11, 2, 11, 0, "Z(MN,EI)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);
  /** Z(MN,EI) Tau(MN,AB) --> W(EI,AB) **/
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WEIAB");
  global_dpd_->contract444(&Z, &T, &W, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);

  /** t_i^f <mn||ef> --> Z(mn,ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 11, 2, 11, 0, "Z(mn,ei)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);
  /** Z(mn,ei) Tau(mn,ab) --> W(ei,ab) **/
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Weiab");
  global_dpd_->contract444(&Z, &T, &W, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);

  /** t_i^f <Mn|Ef> --> Z(Mn,Ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(Mn,Ei)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);
  /** Z(Mn,Ei) Tau(Mn,Ab) --> W(Ei,Ab) **/
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  global_dpd_->contract444(&Z, &T, &W, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);

  /** t_I^F <mN|eF> --> Z(mN,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(mN,eI)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);
  /** Z(mN,eI) Tau(mN,aB) --> W(eI,aB) **/
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
  global_dpd_->contract444(&Z, &T, &W, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);


  /** <MN||EI> Tau(MN,AB) --> W(EI,AB) **/
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 7, 10, 7, 0, "Z(IE,AB)");
  global_dpd_->contract444(&E, &T, &Z, 1, 1, -1.0, 0.0);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 7, "Z(EI,AB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 7, 11, 7, 0, "Z(EI,AB)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WEIAB");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /** <mn||ei> Tau(mn,ab) --> W(ei,ab) **/
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 7, 10, 7, 0, "Z(ie,ab)");
  global_dpd_->contract444(&E, &T, &Z, 1, 1, -1.0, 0.0);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 7, "Z(ei,ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 7, 11, 7, 0, "Z(ei,ab)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Weiab");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /** <Mn|Ei> Tau(Mn,Ab) --> W(Ei,Ab) **/
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(iE,bA)");
  global_dpd_->contract444(&E, &T, &Z, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 5, "Z(Ei,bA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(Ei,bA)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 11, 5, "Z(Ei,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(Ei,Ab)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /** <mN|eI> Tau(mN,aB) --> W(eI,aB) **/
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(Ie,Ba)");
  global_dpd_->contract444(&E, &T, &Z, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 5, "Z(eI,Ba)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(eI,Ba)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 11, 5, "Z(eI,aB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(eI,aB)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /** <MB||EF> t_IM^AF + <MA||FE> t_IM^BF --> W(EI,AB) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, prqs, 10, 5, "F(MF,AE)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "F(MF,AE)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(IB,AE)");
  global_dpd_->contract444(&T2, &F, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 10, 5, "Z(IE,AB)2");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(IE,AB)2");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 5, "Z(EI,AB)2");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(EI,AB)2");
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, pqsr, 11, 5, "Z(EI,BA)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(EI,BA)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "WEIAB");
  global_dpd_->buf4_axpy(&Z1, &W, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&W);

  /** -<mB|fE> t_Im^Af + <mA|fE> t_Im^Bf --> W(EI,AB) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, prqs, 10, 5, "F(mf,AE)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "F(mf,AE)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(IB,AE)");
  global_dpd_->contract444(&T2, &F, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 10, 5, "Z(IE,AB)2");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(IE,AB)2");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 5, "Z(EI,AB)2");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(EI,AB)2");
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, pqsr, 11, 5, "Z(EI,BA)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(EI,BA)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "WEIAB");
  global_dpd_->buf4_axpy(&Z1, &W, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&W);

  /** <mb||ef> t_im^af + <ma||fe> t_im^bf --> W(ei,ab) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, prqs, 10, 5, "F(mf,ae)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "F(mf,ae)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(ib,ae)");
  global_dpd_->contract444(&T2, &F, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 10, 5, "Z(ie,ab)2");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(ie,ab)2");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 5, "Z(ei,ab)2");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(ei,ab)2");
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, pqsr, 11, 5, "Z(ei,ba)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(ei,ba)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "Weiab");
  global_dpd_->buf4_axpy(&Z1, &W, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&W);

  /** -<Mb|Fe> t_iM^aF + <Ma|Fe> t_iM^bF --> W(ei,ab) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, prqs, 10, 5, "F(MF,ae)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "F(MF,ae)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(ib,ae)");
  global_dpd_->contract444(&T2, &F, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, psrq, 10, 5, "Z(ie,ab)2");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(ie,ab)2");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qprs, 11, 5, "Z(ei,ab)2");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(ei,ab)2");
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, pqsr, 11, 5, "Z(ei,ba)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(ei,ba)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "Weiab");
  global_dpd_->buf4_axpy(&Z1, &W, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&W);

  /** -<Mb|Ef> t_Mi^Af - <MA||EF> t_iM^bF + <mA|fE> t_im^bf --> W(Ei,Ab) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, psrq, 10, 5, "F(Mf,Eb)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "F(Mf,Eb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(iA,Eb)");
  global_dpd_->contract444(&T2, &F, &Z, 1, 1, -1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, prqs, 10, 5, "Z(iE,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP1, prqs, 10, 5, "F(mf,AE)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "F(mf,AE)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 10, 5, 10, 5, 0, "Z(ib,AE)");
  global_dpd_->contract444(&T2, &F, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 5, "Z(iE,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP2, psrq, 10, 5, "F(MF,EA)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP2, 0, 10, 5, 10, 5, 0, "F(MF,EA)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 5, 10, 5, 0, "Z(ib,EA)");
  global_dpd_->contract444(&T2, &F, &Z, 0, 1, -1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP4, prqs, 10, 5, "Z(iE,bA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 10, 5, 10, 5, 0, "Z(iE,bA)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, pqsr, 10, 5, "Z(iE,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(iE,Ab)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(iE,Ab)");
  global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 10, 5, 10, 5, 0, "Z(iE,Ab)");
  global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, qprs, 11, 5, "Z(Ei,Ab)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(Ei,Ab)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** -<mB|eF> t_mI^aF - <ma||ef> t_Im^Bf + <Ma|Fe> t_IM^BF --> W(eI,aB) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, psrq, 10, 5, "F(mF,eB)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "F(mF,eB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(Ia,eB)");
  global_dpd_->contract444(&T2, &F, &Z, 1, 1, -1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, prqs, 10, 5, "Z(Ie,aB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP1, prqs, 10, 5, "F(MF,ae)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "F(MF,ae)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 10, 5, 10, 5, 0, "Z(IB,ae)");
  global_dpd_->contract444(&T2, &F, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 5, "Z(Ie,aB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_TMP2, psrq, 10, 5, "F(mf,ea)");
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_init(&F, PSIF_CC_TMP2, 0, 10, 5, 10, 5, 0, "F(mf,ea)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 5, 10, 5, 0, "Z(IB,ea)");
  global_dpd_->contract444(&T2, &F, &Z, 0, 1, -1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP4, prqs, 10, 5, "Z(Ie,Ba)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP4, 0, 10, 5, 10, 5, 0, "Z(Ie,Ba)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, pqsr, 10, 5, "Z(Ie,aB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z(Ie,aB)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 10, 5, 10, 5, 0, "Z(Ie,aB)");
  global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 10, 5, 10, 5, 0, "Z(Ie,aB)");
  global_dpd_->buf4_axpy(&Z1, &Z2, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, qprs, 11, 5, "Z(eI,aB)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(eI,aB)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Final terms of Wabei **/

  /** t_IN^BF <MN||EF> --> Z_IBME **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IB,ME)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->contract444(&T2, &D, &Z, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  global_dpd_->contract444(&T2, &D, &Z, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 10, "Z(IE,MB)");
  global_dpd_->buf4_close(&Z);

  /** t_M^A ( -<MB||IE> + Z1_IEMB ) --> Z2_EIAB **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(IE,MB)");
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  global_dpd_->buf4_axpy(&C, &Z1, -1.0);
  global_dpd_->buf4_close(&C);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z2(IE,AB)");
  global_dpd_->contract244(&T1, &Z1, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, qprs, 11, 5, "Z(EI,AB)2");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(EI,AB)2");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 11, 5, "Z(EI,BA)");
  global_dpd_->buf4_close(&Z);

  /** Z1_EIAB - Z2_EIBA --> W_EIAB **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(EI,AB)2");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(EI,BA)");
  global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "WEIAB");
  global_dpd_->buf4_axpy(&Z2, &W, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&W);


  /** t_in^bf <mn||ef> --> Z_ibme **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ib,me)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->contract444(&T2, &D, &Z, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  global_dpd_->contract444(&T2, &D, &Z, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 10, "Z(ie,mb)");
  global_dpd_->buf4_close(&Z);

  /** t_m^a ( -<mb||ie> + Z1_iemb ) --> Z2_eiab **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(ie,mb)");
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  global_dpd_->buf4_axpy(&C, &Z1, -1.0);
  global_dpd_->buf4_close(&C);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "Z2(ie,ab)");
  global_dpd_->contract244(&T1, &Z1, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, qprs, 11, 5, "Z(ei,ab)2");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(ei,ab)2");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 11, 5, "Z(ei,ba)");
  global_dpd_->buf4_close(&Z);

  /** - Z1_eiab + Z2_eiba --> W_eiab **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z(ei,ab)2");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z(ei,ba)");
  global_dpd_->buf4_axpy(&Z1, &Z2, -1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "Weiab");
  global_dpd_->buf4_axpy(&Z2, &W, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&W);

  /** t_in^bf  <Mn|Ef> + t_iN^bF <MN||EF> --> Z1_MEib **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&Z);

  /** - t_Ni^Af <mN|fE> --> Z2_mEiA **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_sort(&D, PSIF_CC_TMP1, psrq, 10, 11, "D(mE,fN)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "D(mE,fN)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  global_dpd_->buf4_sort(&T2, PSIF_CC_TMP2, psrq, 10, 11, "T2(iA,fN)");
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_init(&D, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "D(mE,fN)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TMP2, 0, 10, 11, 10, 11, 0, "T2(iA,fN)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(mE,iA)");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,ib)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 10, "Z(Mb,iE)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(Mb,iE)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 10, 11, "Z(Mb,Ei)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(mE,iA)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 10, "Z(mA,iE)");
  global_dpd_->buf4_close(&Z);

  /** - t_M^A ( <Mb|Ei> + Z(Mb,Ei) ) --> Z1(Ab,Ei) **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_sort(&D, PSIF_CC_TMP2, psrq, 10, 11, "D(Mb,Ei)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_TMP2, 0, 10, 11, 10, 11, 0, "D(Mb,Ei)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(Mb,Ei)");
  global_dpd_->buf4_axpy(&D, &Z, 1.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z1(Ab,Ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 11, 5, "Z1(Ei,Ab)");
  global_dpd_->buf4_close(&Z1);

  /** t_m^b ( - <mA|iE> + Z(mA,iE) ) --> Z2(Ab,Ei) **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(mA,iE)");
  global_dpd_->buf4_axpy(&C, &Z, -1.0);
  global_dpd_->buf4_close(&C);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, 0, 5, 10, 5, 10, 0, "Z2(bA,iE)");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, qprs, 5, 10, "Z2(Ab,iE)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 5, 10, 5, 10, 0, "Z2(Ab,iE)");
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP2, pqsr, 5, 11, "Z2(Ab,Ei)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z2(Ab,Ei)");
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, rspq, 11, 5, "Z2(Ei,Ab)");
  global_dpd_->buf4_close(&Z2);

  /** Z1(Ei,Ab) + Z2(Ei,Ab) --> W(Ei,Ab) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z1(Ei,Ab)");
  global_dpd_->buf4_axpy(&Z1, &W, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z2(Ei,Ab)");
  global_dpd_->buf4_axpy(&Z2, &W, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&W);

  /** t_IN^BF  <mN|eF> + t_In^Bf <mn||ef> --> Z1_meIB **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(me,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&Z);

  /** - t_Ni^Af <mN|fE> --> Z2_mEiA **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_sort(&D, PSIF_CC_TMP1, psrq, 10, 11, "D(Me,Fn)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "D(Me,Fn)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->buf4_sort(&T2, PSIF_CC_TMP2, psrq, 10, 11, "T2(Ia,Fn)");
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_init(&D, PSIF_CC_TMP1, 0, 10, 11, 10, 11, 0, "D(Me,Fn)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TMP2, 0, 10, 11, 10, 11, 0, "T2(Ia,Fn)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(Me,Ia)");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(me,IB)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 10, "Z(mB,Ie)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(mB,Ie)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 10, 11, "Z(mB,eI)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(Me,Ia)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, psrq, 10, 10, "Z(Ma,Ie)");
  global_dpd_->buf4_close(&Z);

  /** - t_m^a ( <mB|eI> + Z(mB,eI) ) --> Z1(aB,eI) **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_sort(&D, PSIF_CC_TMP2, psrq, 10, 11, "D(mB,eI)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_TMP2, 0, 10, 11, 10, 11, 0, "D(mB,eI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(mB,eI)");
  global_dpd_->buf4_axpy(&D, &Z, 1.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z1(aB,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, rspq, 11, 5, "Z1(eI,aB)");
  global_dpd_->buf4_close(&Z1);

  /** t_M^B ( - <Ma|Ie> + Z(Ma,Ie) ) --> Z2(aB,eI) **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(Ma,Ie)");
  global_dpd_->buf4_axpy(&C, &Z, -1.0);
  global_dpd_->buf4_close(&C);
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, 0, 5, 10, 5, 10, 0, "Z2(Ba,Ie)");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, qprs, 5, 10, "Z2(aB,Ie)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 5, 10, 5, 10, 0, "Z2(aB,Ie)");
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP2, pqsr, 5, 11, "Z2(aB,eI)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, 0, 5, 11, 5, 11, 0, "Z2(aB,eI)");
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP1, rspq, 11, 5, "Z2(eI,aB)");
  global_dpd_->buf4_close(&Z2);

  /** Z1(eI,aB) + Z2(eI,aB) --> W(eI,aB) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "Z1(eI,aB)");
  global_dpd_->buf4_axpy(&Z1, &W, 1.0);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP1, 0, 11, 5, 11, 5, 0, "Z2(eI,aB)");
  global_dpd_->buf4_axpy(&Z2, &W, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&W);
}

}} // namespace psi::cchbar
