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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void V_build(void)
{
  dpdbuf4 V, L, T;
  int G_irr;
  G_irr = params.G_irr;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 2, 2, 2, 2, 0, "VMNIJ");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 2, 2, 2, 2, 0, "Vmnij");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 0, 0, 0, 0, 0, "VMnIj");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 10, 10, 10, 10, 0, "VIAJB");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "LIAJB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 10, 10, 10, 10, 0, "Viajb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "Liajb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 10, 10, 10, 10, 0, "VIAjb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "Liajb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 10, 10, 10, 10, 0, "ViaJB");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "LIAJB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 10, 10, 10, 10, 0, "ViAjB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "LIbjA");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 10, 10, 10, 10, 0, "VIaJb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 10, 10, 10, 10, 0, "LjAIb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 2, 2, 2, 2, 0, "VMNIJ");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 12, 12, 12, 12, 0, "Vmnij");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 12, 17, 12, 17, 0, "Lijab");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 22, 22, 22, 22, 0, "VMnIj");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 20, 20, 20, 20, 0, "VIAJB");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 20, 20, 20, 20, 0, "LIAJB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 30, 30, 30, 30, 0, "Viajb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 30, 30, 30, 30, 0, "Liajb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 20, 30, 20, 30, 0, "VIAjb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 30, 30, 30, 30, 0, "Liajb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 30, 20, 30, 20, 0, "ViaJB");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 20, 20, 20, 20, 0, "LIAJB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 27, 27, 27, 27, 0, "ViAjB");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 24, 27, 24, 27, 0, "LIbjA");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_CC_MISC, G_irr, 24, 24, 24, 24, 0, "VIaJb");
    global_dpd_->buf4_init(&L, PSIF_CC_GLG, G_irr, 27, 24, 27, 24, 0, "LjAIb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);
  }
}

}} // namespace psi::ccdensity
