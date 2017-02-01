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

void CCEnergyWavefunction::WmnijT2(void)
{
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 WMNIJ, Wmnij, WMnIj;
  dpdbuf4 tauIJAB, tauijab, tauIjAb;

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->contract444(&WMNIJ, &tauIJAB, &newtIJAB, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&newtIJAB);

    global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->contract444(&Wmnij, &tauijab, &newtijab, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&newtijab);

    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params_.ref == 2) { /*** UHF ***/

    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->contract444(&WMNIJ, &tauIJAB, &newtIJAB, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&newtIJAB);

    global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->contract444(&Wmnij, &tauijab, &newtijab, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&newtijab);

    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->contract444(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&newtIjAb);

  }
}
}} // namespace psi::ccenergy
