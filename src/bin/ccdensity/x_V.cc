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

/* compute intermediates for excited state pdm that look like
   ground state intermediates, only use GL */

void V_build_x(void)
{
  int L_irr;
  dpdbuf4 V, L, T;
  L_irr = params.L_irr;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAJB");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "Viajb");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_->contract444(&T, &L, &V, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAjb");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 1.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViaJB");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViAjB");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIbjA");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIaJb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LjAIb");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&V);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 20, 20, 20, 20, 0, "VIAJB");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 30, 30, 30, 30, 0, "Viajb");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_->contract444(&T, &L, &V, 1, 1, 1.0, 1.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 20, 30, 20, 30, 0, "VIAjb");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 1.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 30, 20, 30, 20, 0, "ViaJB");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 27, 27, 27, 27, 0, "ViAjB");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 24, 27, 24, 27, 0, "LIbjA");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
    dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&V);

    dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 24, 24, 24, 24, 0, "VIaJb");
    dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 27, 24, 27, 24, 0, "LjAIb");
    dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    dpd_->buf4_close(&T);
    dpd_->buf4_close(&L);
    dpd_->buf4_close(&V);
  }
}

}} // namespace psi::ccdensity
