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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/* WmbejT2(): Contributions of Wmbej intermediates to T2.
**
** t(ij,ab) <--- P(ij) P(ab) t(im,ae) Wmbej
**
** Spin cases for UHF or ROHF orbitals:
** ------------------------------------
**                 *** AA ***
** + T(IA,ME) W(ME,JB) + T(IA,me) W(me,JB)
** - T(JA,ME) W(ME,IB) - T(JA,me) W(me,IB)
** - T(IB,ME) W(ME,JA) - T(IB,me) W(me,JA)
** + T(JB,ME) W(ME,IA) + T(JB,me) W(me,IA)
**
**                 *** BB ***
** + T(ia,me) W(me,jb) + T(ia,ME) W(ME,jb)
** - T(ja,me) W(me,ib) - T(ja,ME) W(ME,ib)
** - T(ib,me) W(me,ja) - T(ib,ME) W(ME,ja)
** + T(jb,me) W(me,ia) + T(jb,ME) W(ME,ia)
**
**                 *** AB ***
** + T(IA,ME) W(ME,jb) + T(IA,me) W(me,jb)
** + T(MA,je) W(Me,Ib) + T(IE,mb) W(mE,jA)
** + T(jb,ME) W(ME,IA) + T(jb,me) W(me,IA)
**
** For the AA and BB spin cases, only the first two terms of each need to
** be evaluated, while for AB, all six terms are different.
**
** The current version of this code requires ten contractions (two each for
** AA and BB and six for AB), one complex sort each for AA and BB, two
** complex sorts for AB, and three simple sorts each for AA and BB.
**
** For RHF orbitals:
** -----------------
** For RHF orbitals, we have two convenient identities:
**
** T(IJ,AB) = T(ij,ab) = T(Ij,Ab) - T(Ij,Ba)
**
** and
**
** W(MB,EJ) = W(mb,ej) = W(Mb,Ej) + W(Mb,eJ)
**
** so that only the AB T2's and the ABAB and ABBA W intermediates are
** necessary.
**
** Therefore, only the AB spin case from above is required for RHF
** orbtials.
**
** + T(IA,ME) W(ME,jb) + T(IA,me) W(me,jb)              I
** + T(MA,je) W(Me,Ib) + T(IE,mb) W(mE,jA)         II   +   III
** + T(jb,ME) W(ME,IA) + T(jb,me) W(me,IA)              IV
**
** The AB T2 term labelled I above may be rewritten as:
**
** 1/2 [ (2 T(IA,me) - T(Ia,mE)) (2 W(ME,jb) + W(Me,Jb))] + 1/2 T(Ia,mE) W(Me,Jb)
**
** Term III from the AB case above is actually the same as the last
** term in the expression above, apart from the 1/2 and swapping of a
** and b indices.  So, for RHF orbitals, we need only evaluate the two
** contractions:
**
**   X(IA,jb) = 1/2 [ (2 T(IA,me) - T(Ia,mE)) (2 W(ME,jb) + W(Me,Jb))]
**
** and
**
**   Z(Ia,Jb) = T(Ia,mE) W(Me,Jb)
**
** Then, I = X(IA,jb) + 1/2 Z(Ia,Jb)
** and   III = Z(Ib,Ja) (i.e., we have to sort Z to get III)
**
** Finally, we obtain II + IV by summing I and III and swapping I<-->j
** and A<-->b.
**
** TDC
** May 2000
** Revised August 2001
** Last revised October 2001
*/

void CCEnergyWavefunction::WmbejT2(void)
{
  dpdbuf4 T2new, T2, W, T2B, W1, W2, Z;

  if(params_.ref == 0) { /** RHF **/
    /*** AB ***/

    /* 2 W(ME,jb) + W(Me,Jb) */
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_copy(&W, PSIF_CC_TMP0, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_init(&W2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->buf4_axpy(&W2, &W1, 2);
    global_dpd_->buf4_close(&W2);
    global_dpd_->buf4_close(&W1);


    /* T2(Ib,mE) * W(mE,jA) --> Z(Ib,jA) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (Ib,jA)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* T2(Ib,jA) --> T2(IA,jb) (part III) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, psrq, 10, 10, "T2 (IA,jb) 3");
    global_dpd_->buf4_close(&T2new);

    /* 1/2 [ (2 T2(IA,me) - T2(IE,ma)) * (2 W(ME,jb) + W(Me,Jb)] --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb) 1");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 0.5, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);
    /* 1/2 Z(Ib,jA) + T2(IA,jb) --> T2(IA,jb) (Part I) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (Ib,jA)");
    global_dpd_->buf4_axpy(&Z, &T2new, 0.5);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&T2new);

    /* T2(IA,jb) (I) + T2(IA,jb) (III) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb) 1");
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "T2 (Ij,Ab) (1+3)");
    global_dpd_->buf4_close(&T2new);

    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) (1+3)");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "T2 (Ij,Ab) (2+4)");
    global_dpd_->buf4_close(&T2new);

    /* T2(Ij,Ab) <--- I + II + III + IV */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) (1+3)");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) (2+4)");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_close(&T2new);

  }
  else if(params_.ref == 1) { /** ROHF **/

    /*** AA ***/

    /* T2(IA,ME) * W(ME,JB) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,JB)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(IA,me) * W(me,JB) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(IA,JB) --> T2(IJ,AB) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);

    /* P(IJ) P(AB) T2(IA,JB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");

    /* T2(IA,JB) - T2(JA,IB) - T2(IB,JA) + T2(JB,IA) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, +1);
    global_dpd_->buf4_close(&T2);

    /* T2(IJ,AB) --> New T2(IJ,AB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /*** BB ***/

    /* T2(ia,me) * W(me,jb) --> T2(ia,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(ia,ME) * W(ME,jb) --> T2(ia,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(ia,jb) --> T2(ij,ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);

    /* P(ij) P(ab) T2(ia,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");

    /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) + T2(ji,ba) --> T2(ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, +1);
    global_dpd_->buf4_close(&T2);

    /* T2(ij,ab) --> New T2(ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tijab");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /*** AB ***/

    /* T2(IA,ME) * W(ME,jb) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (IA,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(IA,me) * W(me,jb) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* W(ME,IA) * T2(jb,ME) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* W(me,IA) * T2(jb,me) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(IA,jb) --> T2(Ij,Ab) (part 1) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "T2 (Ij,Ab) 1");
    global_dpd_->buf4_close(&T2new);

    /* T2(Ib,mE) * W(mE,jA) --> T2(Ib,jA) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "T2 (Ib,jA)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* W(Me,Ib) * T2(jA,Me) --> T2(Ib,jA) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(Ib,jA) --> T2(Ij,Ab) (part 2) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prsq, 0, 5, "T2 (Ij,Ab) 2");
    global_dpd_->buf4_close(&T2new);


    /* T2(Ij,Ab) (part 1) + T2(Ij,Ab) (part 2) --> New T2(Ij,Ab) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) 1");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_close(&T2new);

  } /*** ROHF ***/
  else if(params_.ref == 2) { /*** UHF ***/

    /*** AA ***/

    /* T2(IA,ME) * W(ME,JB) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "T2 (IA,JB)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(IA,me) * W(me,JB) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(IA,JB) --> T2(IJ,AB) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);

    /* P(IJ) P(AB) T2(IA,JB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");

    /* T2(IA,JB) - T2(JA,IB) - T2(IB,JA) + T2(JB,IA) --> T2(IA,JB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, +1);
    global_dpd_->buf4_close(&T2);

    /* T2(IJ,AB) --> New T2(IJ,AB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /*** BB ***/

    /* T2(ia,me) * W(me,jb) --> T2(ia,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "T2 (ia,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(ia,ME) * W(ME,jb) --> T2(ia,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(ia,jb) --> T2(ij,ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 10, 15, "X(10,15) 1");
    global_dpd_->buf4_close(&T2new);

    /* P(ij) P(ab) T2(ia,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 10, 15, "X(10,15) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 10, 15, "X(10,15) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 10, 15, "X(10,15) 4");

    /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) + T2(ji,ba) --> T2(ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, +1);
    global_dpd_->buf4_close(&T2);

    /* T2(ij,ab) --> New T2(ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "New tijab");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /*** AB ***/

    /* T2(IA,ME) * W(ME,jb) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 20, 30, 20, 30, 0, "T2 (IA,jb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(IA,me) * W(me,jb) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* W(ME,IA) * T2(jb,ME) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* W(me,IA) * T2(jb,me) --> T2(IA,jb) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(IA,jb) --> T2(Ij,Ab) (part 1) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prqs, 22, 28, "T2 (Ij,Ab) 1");
    global_dpd_->buf4_close(&T2new);

    /* T2(Ib,mE) * W(mE,jA) --> T2(Ib,jA) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 24, 27, 24, 27, 0, "T2 (Ib,jA)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");
    global_dpd_->contract444(&T2, &W, &T2new, 0, 1, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* W(Me,Ib) * T2(jA,Me) --> T2(Ib,jA) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");
    global_dpd_->contract444(&W, &T2, &T2new, 1, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&T2);

    /* T2(Ib,jA) --> T2(Ij,Ab) (part 2) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, prsq, 22, 28, "T2 (Ij,Ab) 2");
    global_dpd_->buf4_close(&T2new);


    /* T2(Ij,Ab) (part 1) + T2(Ij,Ab) (part 2) --> New T2(Ij,Ab) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "T2 (Ij,Ab) 1");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "T2 (Ij,Ab) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_close(&T2new);

  } /*** UHF ***/
}
}} // namespace psi::ccenergy
