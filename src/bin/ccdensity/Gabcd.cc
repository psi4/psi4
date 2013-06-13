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
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void Gabcd(void)
{
  dpdbuf4 G, L, T;
  int G_irr;
  G_irr = params.G_irr;

  if(params.ref == 0) { /** RHF **/
    dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    if (params.ground)
      dpd_->buf4_symm(&G);
    dpd_->buf4_close(&G);
  }
  else if(params.ref == 1) { /** RHF/ROHF **/

    dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    if (params.ground)
      dpd_->buf4_symm(&G);
    dpd_->buf4_close(&G);

    dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
    dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    if (params.ground)
      dpd_->buf4_symm(&G);
    dpd_->buf4_close(&G);

    dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    if (params.ground)
      dpd_->buf4_symm(&G);
    dpd_->buf4_close(&G);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    dpd_->buf4_symm(&G);
    dpd_->buf4_close(&G);

    dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 17, 17, 17, 17, 0, "Gabcd");
    dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    dpd_->buf4_symm(&G);
    dpd_->buf4_close(&G);

    dpd_->buf4_init(&G, PSIF_CC_GAMMA, G_irr, 28, 28, 28, 28, 0, "GAbCd");
    dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_->contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    dpd_->buf4_symm(&G);
    dpd_->buf4_close(&G);
  }
}

}} // namespace psi::ccdensity
