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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
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
    dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_->contract444(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1.0, 1.0);
    dpd_->buf4_close(&WMnIj);
    dpd_->buf4_close(&LIjAb);
    dpd_->buf4_close(&newLIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_->contract444(&WMNIJ, &LIJAB, &newLIJAB, 0, 1, 1.0, 1.0);
    dpd_->buf4_close(&WMNIJ);
    dpd_->buf4_close(&LIJAB);
    dpd_->buf4_close(&newLIJAB);

    dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    dpd_->contract444(&Wmnij, &Lijab, &newLijab, 0, 1, 1.0, 1.0);
    dpd_->buf4_close(&Wmnij);
    dpd_->buf4_close(&Lijab);
    dpd_->buf4_close(&newLijab);

    dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_->contract444(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1.0, 1.0);
    dpd_->buf4_close(&WMnIj);
    dpd_->buf4_close(&LIjAb);
    dpd_->buf4_close(&newLIjAb);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_->contract444(&WMNIJ, &LIJAB, &newLIJAB, 0, 1, 1, 1);
    dpd_->buf4_close(&WMNIJ);
    dpd_->buf4_close(&LIJAB);
    dpd_->buf4_close(&newLIJAB);

    dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    dpd_->contract444(&Wmnij, &Lijab, &newLijab, 0, 1, 1, 1);
    dpd_->buf4_close(&Wmnij);
    dpd_->buf4_close(&Lijab);
    dpd_->buf4_close(&newLijab);

    dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    dpd_->contract444(&WMnIj, &LIjAb, &newLIjAb, 0, 1, 1, 1);
    dpd_->buf4_close(&WMnIj);
    dpd_->buf4_close(&LIjAb);
    dpd_->buf4_close(&newLIjAb);
  }
}


}} // namespace psi::cclambda
