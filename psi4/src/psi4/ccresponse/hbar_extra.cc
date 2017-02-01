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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void hbar_extra(void) {
  dpdfile2 lt;
  dpdbuf4 W1, W2, W;
  dpdbuf4 t2, l2;

  /* LIjAb * TIjAb */
  global_dpd_->file2_init(&lt, PSIF_CC_OEI, 0, 0, 0, "Lt_IJ");

  global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB 0 -1");
  global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  global_dpd_->contract442(&l2, &t2, &lt, 0, 0, -1.0, 0.0);
  global_dpd_->buf4_close(&t2);
  global_dpd_->buf4_close(&l2);

  global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb 0 -1");
  global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->contract442(&l2, &t2, &lt, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&t2);
  global_dpd_->buf4_close(&l2);

  global_dpd_->file2_close(&lt);

  /* 2 W(ME,jb) + W(Me,Jb) */
  global_dpd_->buf4_init(&W1, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  global_dpd_->buf4_copy(&W1, PSIF_CC_HBAR, "2 W(ME,jb) + W(Me,Jb)");
  global_dpd_->buf4_close(&W1);
  global_dpd_->buf4_init(&W1, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
  global_dpd_->buf4_init(&W2, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  global_dpd_->buf4_axpy(&W2, &W1, 2);
  global_dpd_->buf4_close(&W2);
  global_dpd_->buf4_sort(&W1, PSIF_CC_HBAR, rspq, 10, 10, "2 W(jb,ME) + W(Jb,Me)");
  global_dpd_->buf4_close(&W1);
}

}} // namespace psi::ccresponse
