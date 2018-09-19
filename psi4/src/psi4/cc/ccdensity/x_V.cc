/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
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

/* compute intermediates for excited state pdm that look like
   ground state intermediates, only use GL */

void V_build_x(void)
{
  int L_irr;
  dpdbuf4 V, L, T;
  L_irr = params.L_irr;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAJB");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "Viajb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAjb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViaJB");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViAjB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIbjA");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIaJb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LjAIb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 20, 20, 20, 20, 0, "VIAJB");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 30, 30, 30, 30, 0, "Viajb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 20, 30, 20, 30, 0, "VIAjb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 30, 20, 30, 20, 0, "ViaJB");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->contract444(&T, &L, &V, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 27, 27, 27, 27, 0, "ViAjB");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 24, 27, 24, 27, 0, "LIbjA");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);

    global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 24, 24, 24, 24, 0, "VIaJb");
    global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 27, 24, 27, 24, 0, "LjAIb");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    global_dpd_->contract444(&T, &L, &V, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&V);
  }
}

}} // namespace psi::ccdensity
