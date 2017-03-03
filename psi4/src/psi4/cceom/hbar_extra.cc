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
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void hbar_extra(void) {
  dpdbuf4 W, W1, W2, WAmEf, WmBeJ, WmBEj, WmNIe, WMnIe;

  if (params.eom_ref == 2) {
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 20, 20, 20, 20, 0, "WMBEJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rspq, 20, 20, "WMBEJ (JB,ME)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 30, 20, 30, 20, 0, "WmBeJ"); /* (me,JB) */
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rspq, 20, 30, "WmBeJ (JB,me)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 30, 30, 30, 30, 0, "Wmbej");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rspq, 30, 30, "Wmbej (jb,me)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 20, 30, 20, 30, 0, "WMbEj"); /* (ME,jb) */
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rspq, 30, 20, "WMbEj (jb,ME)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 27, 23, 27, 23, 0, "WmBiJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 27, 22, "WmBiJ (mB,Ji)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 27, 22, 27, 22, 0, "WmBiJ (mB,Ji)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 26, 22, "WmBiJ (Bm,Ji)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 25, 29, 25, 29, 0, "WeIaB");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 24, 29, "WeIaB (Ie,aB)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 24, 29, 24, 29, 0, "WeIaB (Ie,aB)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 24, 28, "WeIaB (Ie,Ab)");
    global_dpd_->buf4_close(&W);
  }

  if(params.eom_ref == 1) {

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMBEJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rspq, 10, 10, "WMBEJ (JB,ME)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WmBeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rspq, 10, 10, "WmBeJ (JB,me)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "Wmbej");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rspq, 10, 10, "Wmbej (jb,me)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, rspq, 10, 10, "WMbEj (jb,ME)");
    global_dpd_->buf4_close(&W);

  }

  if (params.eom_ref == 1) {  /* ROHF */

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WmBiJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 10, 0, "WmBiJ (mB,Ji)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WmBiJ (mB,Ji)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 11, 0, "WmBiJ (Bm,Ji)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WeIaB");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 10, 5, "WeIaB (Ie,aB)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WeIaB (Ie,aB)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 10, 5, "WeIaB (Ie,Ab)");
    global_dpd_->buf4_close(&W);
  }

  if (params.eom_ref == 0 ) { /* RHF */
    /* 2 W(ME,jb) + W(Me,Jb) */
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_copy(&W, PSIF_CC_HBAR, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W1, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_init(&W2, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbEj");
    global_dpd_->buf4_axpy(&W2, &W1, 2);
    global_dpd_->buf4_close(&W2);
    global_dpd_->buf4_sort(&W1, PSIF_CC_HBAR, rspq, 10, 10, "2 W(jb,ME) + W(Jb,Me)");
    global_dpd_->buf4_close(&W1);

    /* used in WamefSD */
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->buf4_scmcopy(&W, PSIF_CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
    global_dpd_->buf4_close(&W);
  }

  return;
}

}} // namespace psi::cceom
