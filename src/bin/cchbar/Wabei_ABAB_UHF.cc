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

/* WAbEi_UHF(): Computes all contributions to the AbEi spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (Ei,Ab) ordering and is referred to on disk as "WAbEi".
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
** For the AbEi spin case, we evaluate these contractions with two
** target orderings, (Ab,Ei) and (Ei,Ab), depending on the term.
** After all terms have been evaluated, the (Ab,Ei) terms are sorted
** into (Ei,Ab) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void WAbEi_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /**** Term I ****/

  /** W(Ei,Ab) <--- <Ei|Ab> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WEiAb");
  global_dpd_->buf4_close(&F);

  /**** Term II ****/

  /** W(Ei,Ab) <--- - F_ME t_Mi^Ab **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1, 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /**** Term III ****/

  /** <Ab|Ef> t_i^f **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&B, &T1, &W, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_close(&W);

  /**** Term IV ****/

  /** WAbEi <-- - t_m^b <mA|fE> t_i^f - t_M^A <Mb|Ef> t_i^f
      Evaluate in three steps: 
          (1) Z_mAEi = - <Am|Ef> t_i^f [stored (Am,Ei)]
	  (2) Z_MbEi = <Mb|Ef> t_i^f   [stored (Mb,Ei)]
          (3) WAbEi <-- t_m^b Z_mAEi - t_M^A Z_MbEi
       Store targets in:  W(Ei,Ab)  and   W'(Ab,Ei)  
  **/

  /** Z(Am,Ei) <-- - <Am|Ef> t_i^f **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, -1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** t_m^b Z(Am,Ei) --> W(Ei,Ab) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Z(Mb,Ei) <-- <Mb|Ef> t_i^f **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** - t_M^A Z(Mb,Ei) --> W'(Ab,Ei) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Term V ****/

  /** WAbEi <-- tau_Mn^Ab <Mn|Ef> t_i^f
      Evaluate in two steps:
         (1) Z_MnEi = <Mn|Ef> t_i^f
         (2) WAbEi <-- tau_Mn^Ab Z_MnEi
      Store target in W'(Ab,Ei)
  **/

  /** Z(Mn,Ei) <-- <Mn|Ef> t_i^f **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 22, 26, 22, 26, 0, "Z(Mn,Ei)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);

  /** tau_Mn^Ab Z1(Mn,Ei) --> W'(Ab,Ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 22, 26, 22, 26, 0, "Z(Mn,Ei)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  global_dpd_->contract444(&T2, &Z, &W, 1, 1, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /**** Term VI ****/

  /** tau_Mn^Ab <Mn|Ei> --> Z(Ab,Ei) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
  global_dpd_->contract444(&T2, &E, &W, 1, 1, 1, 1);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);

  /**** Term VII ****/

  /** WAbEi <-- <bM|fE> t_iM^fA - <AM||EF> t_iM^bF - <Am|Ef> t_im^bf
      Evaluate in five steps:
        (1) Sort <bM|fE> to F(bE,Mf) ordering.  (** Note that we assume a sort has already 
            been done for <AM||EF> and <Am|Ef> in the WABEI and Wabei terms. **)
        (2) Z'(bE,iA) = F(bE,Mf) T(iA,Mf)
        (3) Sort Z'(bE,iA) --> W(Ei,Ab)
        (4) Z''(AE,ib) = - F(AE,MF) T(ib,MF) - F(AE,mf) T(ib,mf)
	(5) Sort Z''(AE,ib) --> W(Ei,Ab)

      NB: The storage for the sorts is expensive and will eventually require out-of-core 
          codes.
  **/

  /** <bM|fE> --> F(bE,Mf) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, psqr, 29, 24, "F <aI|bC> (aC,Ib)");
  global_dpd_->buf4_close(&F);

  /** <bM|fE> t_iM^fA --> Z(bE,iA) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 29, 27, 29, 27, 0, "Z(bE,iA)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 29, 24, 29, 24, 0, "F <aI|bC> (aC,Ib)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, -1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z(bE,iA) --> W(Ei,Ab) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 29, 27, 29, 27, 0, "Z(bE,iA)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, qrsp, 26, 28, "WEiAb", 1);
  global_dpd_->buf4_close(&Z);

  /** Z''(AE,ib) <-- - <AM||EF> t_iM^bF **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 30, 5, 30, 0, "Z(AE,ib)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 5, 20, 5, 20, 0, "F <AI||BC> (AB,IC)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z''(AE,ib) <-- -<Am|Ef> t_im^bf **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 30, 5, 30, 0, "Z(AE,ib)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 5, 30, 5, 30, 0, "F <Ai|Bc> (AB,ic)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z''(AE,ib) --> W(Ei,Ab) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 30, 5, 30, 0, "Z(AE,ib)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, qrps, 26, 28, "WEiAb", 1);
  global_dpd_->buf4_close(&Z);

  /**** Terms VIII and IX ****/

  /** WAbEi <-- - t_M^A { <Mb|Ei> + t_iN^bF <MN||EF> + t_in^bf <Mn|Ef> }
                + t_m^b {-<mA|iE> + t_iN^fA <mN|fE> } 
      Evaluate in three steps: 
         (1) Z_MbEi =  <Mb|Ei> + t_iN^bF <MN||EF> + tin^bf <Mn|Ef>  [stored (Mb,Ei)]
         (2) Z_mAEi = -<mA|iE> + t_iN^fA <mN|fE>                    [stored (Am,Ei)]
         (3) WAbEi <-- - t_M^A Z_MbEi + t_m^b Z_mAEi
      Store targets in     W'(Ab,Ei) and  W(Ei,AB)
  **/

  /** Z(Mb,Ei) <-- <Mb|Ei> **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
  global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "Z(Mb,Ei)");
  global_dpd_->buf4_close(&D);

  /** <MN||EF> t_iN^bF --> Z(ME,ib) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <Mn|Ef> t_in^bf --> Z(ME,ib) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(ME,ib) --> Z(Mb,Ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,ib)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, psqr, 24, 26, "Z(Mb,Ei)", 1);
  global_dpd_->buf4_close(&Z);

  /** W'(Ab,Ei) <-- - t_M^A Z(Mb,Ei) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Z(Am,Ei) <-- - <mA|iE> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 26, 26, 26, 26, 0, "C <Ai|Bj>");
  global_dpd_->buf4_copy(&C, PSIF_CC_TMP0, "Z(Am,Ei)");
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  global_dpd_->buf4_scm(&Z, -1);
  global_dpd_->buf4_close(&Z);

  /** Z(mE,iA) <-- t_iN^fA <mN|fE> **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mE,iA)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(mE,iA) --> Z(Am,Ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mE,iA)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, spqr, 26, 26, "Z(Am,Ei)", 1);
  global_dpd_->buf4_close(&Z);

  /** W(Ei,AB) <-- t_m^b Z_mAEi **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Combine accumulated W'(Ab,Ei) and W(Ei,Ab) terms into WEiAb ****/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rspq, 26, 28, "WEiAb", 1);
  global_dpd_->buf4_close(&W);
}

}} // namespace psi::cchbar