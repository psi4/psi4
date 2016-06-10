/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* WaBeI_UHF(): Computes all contributions to the aBeI spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (eI,aB) ordering and is referred to on disk as "WaBeI".
**
** The spin-orbital expression for the Wabei elements is:
**
** Wabei = <ab||ei> - Fme t_mi^ab + t_i^f <ab||ef>
**         - P(ab) t_m^b <am||ef>t_i^f + 1/2 tau_mn^ab <mn||ef> t_i^f
**         + 1/2 <mn||ei> tau_mn^ab - P(ab) <mb||ef> t_mi^af
**         - P(ab) t_m^a { <mb||ei> - t_ni^bf <mn||ef> }
**
** (cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).)
**
** For the aBeI spin case, we evaluate these contractions with two
** target orderings, (aB,eI) and (eI,aB), depending on the term.
** After all terms have been evaluated, the (aB,eI) terms are sorted
** into (eI,aB) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void WaBeI_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /**** Term I ****/

  /** W(eI,aB) <--- <eI|aB> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WeIaB");
  global_dpd_->buf4_close(&F);

  /**** Term II ****/

  /** W(eI,aB) <--- - F_me t_mI^aB **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1, 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /**** Term III ****/  /** This will require special out-of-core code **/

  /** Z(Ie,Ba) <--- t_I^F <Fe|Ba> **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 28, 24, 28, 0, "Z(Ie,Ba)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &B, &Z, 1, 0, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  /** Z(Ie,Ba) --> W'(aB,eI) **/
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, srqp, 29, 25, "W'(aB,eI)");
  global_dpd_->buf4_close(&Z);

  /**** Term IV ****/

  /** WaBeI <-- - t_M^B <Ma|Fe> t_I^F - t_m^a <mB|eF> t_I^F
      Evaluate in three steps: 
          (1) Z_MaeI = - <aM|eF> t_I^F [stored (aM,eI)]
	  (2) Z_mBeI = <mB|eF> t_I^F   [stored (mB,eI)]
          (3) WaBeI <-- t_M^B Z_MaeI - t_m^a Z_mBeI
       Store targets in:  W(eI,aB)  and   W'(aB,eI)  
  **/

  /** Z(aM,eI) <-- - <aM|eF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, -1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** t_M^B Z(aM,eI) --> W(eI,aB) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Z(mB,eI) <-- <mB|eF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** - t_m^a Z(mB,eI) --> W'(aB,eI) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Term V ****/

  /** WaBeI <-- tau_mN^aB <mN|eF> t_I^F
      Evaluate in two steps:
         (1) Z_mNeI = <mN|eF> t_I^F
         (2) WaBeI <-- tau_mN^aB Z_mNeI
      Store target in W'(aB,eI)
  **/

  /** Z(mN,eI) <-- <mN|eF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 23, 25, 23, 25, 0, "Z(mN,eI)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);

  /** tau_mN^aB Z(mN,eI) --> W'(aB,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 23, 25, 23, 25, 0, "Z(mN,eI)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
  global_dpd_->contract444(&T2, &Z, &W, 1, 1, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /**** Term VI ****/

  /** tau_mN^aB <mN|eI> --> W'(aB,eI) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
  global_dpd_->contract444(&T2, &E, &W, 1, 1, 1, 1);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);

  /**** Term VII ****/

  /** WaBeI <-- <Bm|Fe> t_Im^Fa - <am||ef> t_Im^Bf - <aM|eF> t_IM^BF
      Evaluate in five steps:
        (1) Sort <Bm|Fe> to F(Be,mF) ordering.  (** Note that we assume a sort has already 
            been done for <am||ef> and <aM|eF> in the WABEI and Wabei terms. **)
        (2) Z'(Be,Ia) = F(Be,mF) T(Ia,mF)
        (3) Sort Z'(Be,Ia) --> W(eI,aB)
        (4) Z''(ae,IB) = - F(ae,mf) T(IB,mf) - F(ae,MF) T(IB,MF)
	(5) Sort Z''(ae,IB) --> W(eI,aB)

      NB: The storage for the sorts is expensive and will eventually require out-of-core 
          codes.
  **/

  /** <Bm|Fe> --> F(Be,mF) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, psqr, 28, 27, "F <Ai|Bc> (Ac,iB)");
  global_dpd_->buf4_close(&F);

  /** <Bm|Fe> t_Im^Fa --> Z(Be,Ia) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 28, 24, 28, 24, 0, "Z(Be,Ia)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 28, 27, 28, 27, 0, "F <Ai|Bc> (Ac,iB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, -1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z(Be,Ia) --> W(eI,aB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 28, 24, 28, 24, 0, "Z(Be,Ia)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, qrsp, 25, 29, "WeIaB", 1);
  global_dpd_->buf4_close(&Z);

  /** Z''(ae,IB) <-- - <am||ef> t_Im^Bf **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0, "Z(ae,IB)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 30, 15, 30, 0, "F <ai||bc> (ab,ic)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z''(ae,IB) <-- -<aM|eF> t_IM^BF **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0, "Z(ae,IB)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 20, 15, 20, 0, "F <aI|bC> (ab,IC)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z''(ai,IB) --> W(eI,aB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0, "Z(ae,IB)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, qrps, 25, 29, "WeIaB", 1);
  global_dpd_->buf4_close(&Z);


  /**** Terms VIII and IX ****/

  /** WaBeI <-- - t_m^a { <mB|eI> + t_In^Bf <mn||ef> + t_IN^BF <mN|eF> }
                + t_M^B {-<Ma|Ie> + t_In^Fa <Mn|Fe> } 
      Evaluate in three steps: 
         (1) Z_mBeI =  <mB|eI> + t_In^Bf <mn||ef> + tIN^BF <mN|eF>  [stored (mB,eI)]
         (2) Z_MaeI = -<Ma|Ie> + t_In^Fa <Mn|Fe>                    [stored (aM,eI)]
         (3) WaBeI <-- - t_m^a Z_mBeI + t_M^B Z_MaeI
      Store targets in     W'(aB,eI) and  W(eI,aB)
  **/

  /** Z(mB,eI) <-- <mB|eI> **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
  global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "Z(mB,eI)");
  global_dpd_->buf4_close(&D);

  /** <mn||ef> t_In^Bf --> Z(me,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <mN|eF> t_IN^BF --> Z(me,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(me,IB) --> Z(mB,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, psqr, 27, 25, "Z(mB,eI)", 1);
  global_dpd_->buf4_close(&Z);

  /** W'(aB,eI) <-- - t_m^a Z(mB,eI) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Z(aM,eI) <-- - <Ma|Ie> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
  global_dpd_->buf4_sort(&C, PSIF_CC_TMP0, qpsr, 25, 25, "Z(aM,eI)");
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->buf4_scm(&Z, -1);
  global_dpd_->buf4_close(&Z);

  /** Z(Me,Ia) <-- t_In^Fa <Mn|Fe> **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(Me,Ia) --> Z(aM,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, spqr, 25, 25, "Z(aM,eI)", 1);
  global_dpd_->buf4_close(&Z);

  /** W(eI,aB) <-- t_M^B Z_MaeI **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Combine accumulated W'(aB,eI) and W(eI,aB) terms into WeIaB ****/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rspq, 25, 29, "WeIaB", 1);
  global_dpd_->buf4_close(&W);
}

}} // namespace psi::cchbar