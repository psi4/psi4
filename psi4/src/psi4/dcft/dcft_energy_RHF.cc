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

#include "dcft.h"
#include "psi4/psifiles.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "defines.h"

namespace psi{ namespace dcft{

/**
 * Uses the intermediates to compute the energy for RHF reference
 */
void
DCFTSolver::compute_dcft_energy_RHF()
{
    /*
     * E = lambda_IjAb * ( 2 M_IjAb - M_JiAb )
     * where M_IjAb = G_IjAb + gbar_IjAb
     */
    dcft_timer_on("DCFTSolver::compute_dcft_energy()");

    dpdbuf4 L, G, M, temp;
    double E_test;
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>"); // Lambda <Oo|Vv>

    dcft_timer_on("DCFTSolver::G_IjAb + g_IjAb");
    // M_IjAb = G_IjAb
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>"); // G <Oo|Vv>
    global_dpd_->buf4_copy(&G, PSIF_DCFT_DPD, "M <OO|VV>"); // M <Oo|Vv>
    global_dpd_->buf4_close(&G);
    // M_IjAb += gbar_IjAb
    global_dpd_->buf4_init(&M, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "M <OO|VV>"); // M <Oo|Vv>
    global_dpd_->buf4_init(&temp, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>"); // MO Ints <Oo|Vv>
    dpd_buf4_add(&M, &temp, 1.0);
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_close(&temp);
    dcft_timer_off("DCFTSolver::G_IjAb + g_IjAb");

    // Form (2 M_IjAb - M_JiAb) = (M_IjAb - M_JiAb) + M_IjAb
    global_dpd_->buf4_init(&M, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 1, "M <OO|VV>");
    global_dpd_->buf4_copy(&M, PSIF_DCFT_DPD, "M(temp) <OO|VV>");
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_init(&M, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "M(temp) <OO|VV>");
    global_dpd_->buf4_init(&temp, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "M <OO|VV>");
    dpd_buf4_add(&M, &temp, 1.0);

    E_test = global_dpd_->buf4_dot(&L, &M);
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_close(&temp);
    global_dpd_->buf4_close(&L);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    lambda_energy_ = E_test;

    dcft_timer_off("DCFTSolver::compute_dcft_energy()");

}

}}
