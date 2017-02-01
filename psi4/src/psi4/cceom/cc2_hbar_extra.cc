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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void cc2_hbar_extra(void) {
  dpdbuf4 W1, W2, W;

  if(params.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&W1, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Me,Jb)");
    global_dpd_->buf4_copy(&W1, PSIF_CC2_HET1, "CC2 2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_close(&W1);
    global_dpd_->buf4_init(&W1, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_init(&W2, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj (ME,jb)");
    global_dpd_->buf4_axpy(&W2, &W1, 2);
    global_dpd_->buf4_close(&W2);
    global_dpd_->buf4_close(&W1);
    global_dpd_->buf4_init(&W1, PSIF_CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->buf4_sort(&W1, PSIF_CC2_HET1, rspq, 10, 10, "CC2 2 W(jb,ME) + W(Jb,Me)");
    global_dpd_->buf4_close(&W1);

    /*
    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    dpd_buf4_scmcopy(&W, CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
    dpd_buf4_sort_axpy(&W, CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
    dpd_buf4_sort(&W, CC2_HET1, rspq, 11, 5, "CC2 WAbEi (Ei,Ab)");
    dpd_buf4_close(&W);
    */
  }
}

}} // namespace psi::cceom
