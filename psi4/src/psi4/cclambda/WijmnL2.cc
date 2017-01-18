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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void WijmnL2(int L_irr)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdbuf4 WMNIJ, Wmnij, WMnIj;

  /* RHS += Lmnab*Wijmn */
  if(params.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->contract444(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&LIjAb);
    global_dpd_->buf4_close(&newLIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->contract444(&WMNIJ, &LIJAB, &newLIJAB, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&LIJAB);
    global_dpd_->buf4_close(&newLIJAB);

    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    global_dpd_->contract444(&Wmnij, &Lijab, &newLijab, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&Lijab);
    global_dpd_->buf4_close(&newLijab);

    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->contract444(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&LIjAb);
    global_dpd_->buf4_close(&newLIjAb);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->contract444(&WMNIJ, &LIJAB, &newLIJAB, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&LIJAB);
    global_dpd_->buf4_close(&newLIJAB);

    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    global_dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    global_dpd_->contract444(&Wmnij, &Lijab, &newLijab, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&Lijab);
    global_dpd_->buf4_close(&newLijab);

    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    global_dpd_->contract444(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1, 1);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&LIjAb);
    global_dpd_->buf4_close(&newLIjAb);
  }
}


}} // namespace psi::cclambda
