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

/* WABEI_UHF(): Computes all contributions to the ABEI spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (EI,AB) ordering and is referred to on disk as "WEIAB".
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
** For the ABEI spin case, we evaluate these contractions with two
** target orderings, (AB,EI) and (EI,AB), depending on the term.
** After all terms have been evaluated, the (AB,EI) terms are sorted
** into (EI,AB) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void WABEI_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /**** Term I ****/

  /** W(EI,AB) <--- <EI||AB> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WEIAB");
  global_dpd_->buf4_close(&F);

  /**** Term II ****/

  /** W(EI,AB) <--- - F_ME t_MI^AB **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /**** Term III ****/

  /** W'(AB,EI) <--- <AB||EF> t_I^F **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&B, &T1, &W, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_close(&W);

  /**** Term IV ****/

  /** WABEI <-- t_M^B <MA||EF> t_I^F - t_M^A <MB||EF> t_I^F
      Evaluate in two steps: 
          (1) Z_MBEI = <MB||EF> t_I^F 
          (2) WABEI <-- t_M^B Z_MAEI - t_M^A Z_MBEI
  **/

  /** Z(MB,EI) <-- - <MB||EF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** t_M^A Z(MB,EI) --> Z1(AB,EI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, qprs, 5, 21, "Z2(BA,EI)");
  global_dpd_->buf4_close(&Z1);

  /** Z1(AB,EI) - Z2(BA,EI) --> W'(AB,EI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z2(BA,EI)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 5, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_axpy(&Z1, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z1);

  /**** Term V ****/
  
  /** WABEI <-- 1/2 tau_MN^AB <MN||EF> t_I^F
      Evaluate in two steps:
         (1) Z_MNEI = <MN||EF> t_I^F
         (2) WABEI <-- 1/2 tau_MN^AB Z_MNEI
      Store target in W'(AB,EI)
  **/

  /** Z(MN,EI) <-- <MN||EF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 21, 2, 21, 0, "Z(MN,EI)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** tau_MN^AB Z(MN,EI) --> W'(AB,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 21, 2, 21, 0, "Z(MN,EI)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->contract444(&T2, &Z, &W, 1, 1, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /**** Term VI ****/

  /** tau_MN^AB <MN||EI> --> W'(AB,EI) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
  global_dpd_->contract444(&T2, &E, &W, 1, 1, -1, 1);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);

  /**** Term VII ****/

  /** WABEI <-- <BM||EF> t_IM^AF + <Bm|Ef> t_Im^Af - <AM||EF> t_IM^BF - <Am|Ef> t_Im^Bf
      Evaluate in six steps:
        (1) Sort <BM||EF> and <Bm|Ef> to F(BE,MF) and F(BE,mf) ordering.
        (2) Z(BE,IA) = F(BE,MF) T(IA,MF) + F(BE,mf) T(IA,mf)
        (3) Sort Z(BE,IA) --> Z'(EI,AB)
	(4) Sort Z'(EI,AB) --> Z''(EI,BA)
        (5) AXPY: Z'(EI,AB) = Z'(EI,AB) - Z''(EI,BA)
        (6) AXPY: W(EI,AB) <-- Z'(EI,AB)
      NB: The storage for the sorts is expensive and will eventually require out-of-core 
          codes.
  **/

  /** <BM||EF> --> F(BE,MF) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 5, 21, 5, 1, "F <AI|BC>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 5, 20, "F <AI||BC> (AB,IC)");
  global_dpd_->buf4_close(&F);

  /** <Bm|Ef> --> (BE,mf) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 5, 30, "F <Ai|Bc> (AB,ic)");
  global_dpd_->buf4_close(&F);

  /** <BM||EF> t_IM^AF --> Z(BE,IA) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(BE,IA)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 5, 20, 5, 20, 0, "F <AI||BC> (AB,IC)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, -1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** <Bm|Ef> t_Im^Af --> Z(BE,IA) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(BE,IA)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 5, 30, 5, 30, 0, "F <Ai|Bc> (AB,ic)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, -1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z(BE,IA) --> Z'(EI,AB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(BE,IA)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qrsp, 21, 5, "Z'(EI,AB)");
  global_dpd_->buf4_close(&Z);

  /** Z'(EI,AB) --> Z''(EI,BA) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "Z'(EI,AB)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 21, 5, "Z''(EI,BA)");
  global_dpd_->buf4_close(&Z);

  /** Z'(EI,AB) = Z'(EI,AB) - Z''(EI,BA) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "Z'(EI,AB)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "Z''(EI,BA)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Z1);

  /** W(EI,AB) <-- Z'(EI,AB) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 5, 21, 7, 0, "WEIAB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "Z'(EI,AB)");
  global_dpd_->buf4_axpy(&Z, &W, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Terms VIII and IX ****/

  /** WABEI <-- -P(AB) t_M^A { <MB||EI> + t_IN^BF <MN||EF> + t_In^Bf <Mn|Ef> }
      Evaluate in two steps: 
         (1) Z_MBEI = <MB||EI> + t_IN^BF <MN||EF> + tIn^Bf <Mn|Ef>
         (2) WABEI <-- - t_M^A Z_MBEI + t_M^B Z_MAEI
      Store target in W'(AB,EI)
  **/

  /** Z(MB,EI) <-- <MB||EI> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
  global_dpd_->buf4_copy(&C, PSIF_CC_TMP0, "Z(MB,EI)");
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  global_dpd_->buf4_scm(&Z, -1);
  global_dpd_->buf4_close(&Z);

  /** <MN||EF> t_IN^BF --> Z(ME,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <Mn|Ef> t_In^Bf --> Z(ME,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(ME,IB) --> Z(MB,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, psqr, 20, 21, "Z(MB,EI)", 1);
  global_dpd_->buf4_close(&Z);

  /** Z(AB,EI) <-- -t_M^A Z(MB,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z1, &Z, 0, 0, 0, -1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&Z);

  /** Z(AB,EI) --> Z'(BA,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qprs, 5, 21, "Z'(BA,EI)");
  global_dpd_->buf4_close(&Z);

  /** Z(AB,EI) = Z(AB,EI) - Z'(BA,EI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z'(BA,EI)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Z1);

  /** Z(AB,EI) --> W'(AB,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 5, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /**** Combine accumulated W'(AB,EI) and W(EI,AB) terms into WEIAB ****/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rspq, 21, 7, "WEIAB", 1);
  global_dpd_->buf4_close(&W);
}

}} // namespace psi::cchbar
