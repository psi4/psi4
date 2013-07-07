/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
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

/* Wabei_UHF(): Computes all contributions to the abei spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (ei,ab) ordering and is referred to on disk as "Wabei".
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
** For the abei spin case, we evaluate these contractions with two
** target orderings, (ab,ei) and (ei,ab), depending on the term.
** After all terms have been evaluated, the (ab,ei) terms are sorted
** into (ei,ab) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void Wabei_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /**** Term I ****/

  /** W(ei,ab) <--- <ei||ab> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "Weiab");
  global_dpd_->buf4_close(&F);

  /**** Term II ****/

  /** W(ei,ab) <--- - F_me t_mi^ab **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 31, 17, 31, 17, 0, "Weiab");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /**** Term III ****/

  /** <ab||ef> t_i^f **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 17, 31, 17, 31, 0, "W'(ab,ei)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 17, 15, 15, 15, 1, "B <ab|cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&B, &T1, &W, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_close(&W);

  /**** Term IV ****/

  /** Wabei <-- t_m^b <ma||ef> t_i^f - t_m^a <mb||ef> t_i^f
      Evaluate in two steps: 
          (1) Z_mbei = <mb||ef> t_i^f 
          (2) Wabei <-- t_m^b Z_maei - t_m^a Z_mbei
  **/

  /** Z(mb,ei) <-- - <mb||ef> t_i^f **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Z(mb,ei)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** t_m^a Z(mb,ei) --> Z1(ab,ei) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "Z1(ab,ei)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Z(mb,ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, qprs, 15, 31, "Z2(ba,ei)");
  global_dpd_->buf4_close(&Z1);

  /** Z1(ab,ei) - Z2(ba,ei) --> Z(ab,ei) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "Z1(ab,ei)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "Z2(ba,ei)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 15, 31, 17, 31, 0, "W'(ab,ei)");
  global_dpd_->buf4_axpy(&Z1, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z1);

  /**** Term V ****/

  /** Wabei <-- 1/2 tau_mn^ab <mn||ef> t_i^f
      Evaluate in two steps:
         (1) Z_mnei = <mn||ei> t_i^f
         (2) Wabei <-- 1/2 tau_mn^ab Z_mnei
      Store target in W'(ab,ei)
  **/

  /** Z(mn,ei) <-- <mn||ef> t_i^f **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 12, 31, 12, 31, 0, "Z(mn,ei)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** tau_mn^ab Z(mn,ei) --> W'(ab,ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 12, 31, 12, 31, 0, "Z(mn,ei)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 17, 31, 17, 31, 0, "W'(ab,ei)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
  global_dpd_->contract444(&T2, &Z, &W, 1, 1, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /**** Term VI ****/

  /** tau_mn^ab <mn||ei> --> W'(ab,ei) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 17, 31, 17, 31, 0, "W'(ab,ei)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 12, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
  global_dpd_->contract444(&T2, &E, &W, 1, 1, -1, 1);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);

  /**** Term VII ****/

  /** Wabei <-- <bm||ef> t_im^af + <bM|eF> t_iM^aF - <am||ef> t_im^bf - <aM|eF> t_iM^bF
      Evaluate in six steps:
        (1) Sort <bm||ef> and <bM|eF> to F(be,mf) and F(be,MF) ordering.
        (2) Z(be,ia) = F(be,mf) T(ia,mf) + F(be,MF) T(ia,MF)
        (3) Sort Z(be,ia) --> Z'(ei,ab)
	(4) Sort Z'(ei,ab) --> Z''(ei,ba)
        (5) AXPY: Z'(ei,ab) = Z'(ei,ab) - Z''(ei,ba)
        (6) AXPY: W(ei,ab) <-- Z'(ei,ab)
      NB: The storage for the sorts is expensive and will eventually require out-of-core 
          codes.
  **/

  /** <bm||ef> --> F(be,mf) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 15, 31, 15, 1, "F <ai|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 15, 30, "F <ai||bc> (ab,ic)");
  global_dpd_->buf4_close(&F);

  /** <bM|eF> --> (be,MF) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 15, 20, "F <aI|bC> (ab,IC)");
  global_dpd_->buf4_close(&F);

  /** <bm||ef> t_im^af --> Z(be,ia) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 30, 15, 30, 0, "Z(be,ia)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 30, 15, 30, 0, "F <ai||bc> (ab,ic)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, -1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** <bm|eF> t_iM^aF --> Z(be,ia) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 30, 15, 30, 0, "Z(be,ia)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 20, 15, 20, 0, "F <aI|bC> (ab,IC)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, -1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z(be,ia) --> Z'(ei,ab) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 30, 15, 30, 0, "Z(be,ia)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qrsp, 31, 15, "Z'(ei,ab)");
  global_dpd_->buf4_close(&Z);

  /** Z'(ei,ab) --> Z''(ei,ba) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 31, 15, 31, 15, 0, "Z'(ei,ab)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 31, 15, "Z''(ei,ba)");
  global_dpd_->buf4_close(&Z);

  /** Z'(ei,ab) = Z'(ei,ab) - Z''(ei,ba) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 31, 15, 31, 15, 0, "Z'(ei,ab)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 31, 15, 31, 15, 0, "Z''(ei,ba)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Z1);

  /** W(ei,ab) <-- Z'(ei,ab) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 31, 15, 31, 17, 0, "Weiab");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 31, 15, 31, 15, 0, "Z'(ei,ab)");
  global_dpd_->buf4_axpy(&Z, &W, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Terms VIII and IX ****/

  /** Wabei <-- -P(ab) t_m^a { <mb||ei> + t_in^bf <mn||ef> + t_iN^bF <mN|eF> }
      Evaluate in two steps: 
         (1) Z_mbei = <mb||ei> + t_in^bf <mn||ef> + tiN^bF <mN|eF>
         (2) Wabei <-- - t_m^a Z_mbei + t_m^b Z_maei
      Store target in W'(ab,ei)
  **/

  /** Z(mb,ei) <-- <mb||ei> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
  global_dpd_->buf4_copy(&C, PSIF_CC_TMP0, "Z(mb,ei)");
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Z(mb,ei)");
  global_dpd_->buf4_scm(&Z, -1);
  global_dpd_->buf4_close(&Z);

  /** <mn||ef> t_in^bf --> Z(me,ib) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(me,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <mN|eF> t_iN^bF --> Z(me,ib) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(me,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(me,ib) --> Z(mb,ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(me,ib)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, psqr, 30, 31, "Z(mb,ei)", 1);
  global_dpd_->buf4_close(&Z);

  /** Z(ab,ei) <-- -t_m^a Z(mb,ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "Z(ab,ei)");
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 30, 31, 30, 31, 0, "Z(mb,ei)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z1, &Z, 0, 0, 0, -1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&Z);

  /** Z(ab,ei) --> Z'(ba,ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "Z(ab,ei)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qprs, 15, 31, "Z'(ba,ei)");
  global_dpd_->buf4_close(&Z);

  /** Z(ab,ei) = Z(ab,ei) - Z'(ba,ei) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "Z(ab,ei)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "Z'(ba,ei)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Z1);

  /** Z(ab,ei) --> W'(ab,ei) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 31, 15, 31, 0, "Z(ab,ei)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 15, 31, 17, 31, 0, "W'(ab,ei)");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /**** Combine accumulated W'(ab,ei) and W(ei,ab) terms into Weiab ****/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 17, 31, 17, 31, 0, "W'(ab,ei)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rspq, 31, 17, "Weiab", 1);
  global_dpd_->buf4_close(&W);
}

}} // namespace psi::cchbar
