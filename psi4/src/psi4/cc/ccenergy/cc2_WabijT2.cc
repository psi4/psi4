/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::cc2_WabijT2() {
    dpdbuf4 W;

    if (params_.ref == 0) { /*** RHF ***/

        global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 0, 5, 0, 5, 0, "CC2 WAbIj (Ij,Ab)");
        global_dpd_->buf4_copy(&W, PSIF_CC_TAMPS, "New tIjAb");
        global_dpd_->buf4_close(&W);

    }

    else if (params_.ref == 1) { /*** ROHF ***/

        global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 2, 7, 2, 7, 0, "CC2 Wabij (i>j,a>b)");
        global_dpd_->buf4_copy(&W, PSIF_CC_TAMPS, "New tIJAB");
        global_dpd_->buf4_copy(&W, PSIF_CC_TAMPS, "New tijab");
        global_dpd_->buf4_close(&W);

        global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 0, 5, 0, 5, 0, "CC2 WAbIj (Ij,Ab)");
        global_dpd_->buf4_copy(&W, PSIF_CC_TAMPS, "New tIjAb");
        global_dpd_->buf4_close(&W);

    }

    else if (params_.ref == 2) { /*** UHF ***/

        global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 2, 7, 2, 7, 0, "CC2 WABIJ (I>J,A>B)");
        global_dpd_->buf4_copy(&W, PSIF_CC_TAMPS, "New tIJAB");
        global_dpd_->buf4_close(&W);

        global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 12, 17, 12, 17, 0, "CC2 Wabij (i>j,a>b)");
        global_dpd_->buf4_copy(&W, PSIF_CC_TAMPS, "New tijab");
        global_dpd_->buf4_close(&W);

        global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 22, 28, 22, 28, 0, "CC2 WAbIj (Ij,Ab)");
        global_dpd_->buf4_copy(&W, PSIF_CC_TAMPS, "New tIjAb");
        global_dpd_->buf4_close(&W);
    }
}
}  // namespace ccenergy
}  // namespace psi
