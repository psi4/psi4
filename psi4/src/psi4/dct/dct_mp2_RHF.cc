/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include "dct.h"
#include "psi4/psifiles.h"

#include "psi4/psifiles.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

#include <vector>

namespace psi {
namespace dct {
/**
 * Computes the Hartree-Fock energy and then the MP2 energy as an initial guess.
 * This code is responible for initializing the integral transformation too.
 */
void DCTSolver::initialize_amplitudes_RHF() {
    dct_timer_on("DCTSolver::initialize_amplitudes()");

    std::string guess = options_.get_str("DCT_GUESS");

    outfile->Printf("\tComputing MP2 amplitude guess...\n\n");

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 I, D;

    /*
     * In spin-adapted closed-shell system, only alpha-beta case is needed for computing energy
     */

    // L_IjAb = <Ij|Ab> / D_IjAb
    dct_timer_on("DCTSolver::g_IJAB / D_IJAB");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "MO Ints <OO|VV>");                         // MO Ints <Oo|Vv>
    global_dpd_->buf4_copy(&I, PSIF_DCT_DPD, "Amplitude SF <OO|VV>");  // Amplitude <Oo|Vv>
    global_dpd_->buf4_close(&I);

    global_dpd_->buf4_init(&D, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0,
                           "D <OO|VV>");  // D <Oo|Vv>
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");  // Amplitude <Oo|Vv>
    global_dpd_->buf4_dirprd(&D, &I);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&D);
    dct_timer_off("DCTSolver::g_IJAB / D_IJAB");

    /* build lambda <OO|VV> for tau and G intermediates */
    dpdbuf4 T;
    // Amplitude_IJAB = Amplitude_IjAb - Amplitude_JiAb
    global_dpd_->buf4_init(&I, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1,
                           "Amplitude SF <OO|VV>");
    global_dpd_->buf4_copy(&I, PSIF_DCT_DPD, "Amplitude <OO|VV>");
    // The purpose of having Amplitude <oo|vv> is for better performance of DIIS
    global_dpd_->buf4_copy(&I, PSIF_DCT_DPD, "Amplitude <oo|vv>");
    global_dpd_->buf4_close(&I);

    /*
     * E = lambda_IjAb * M_IjAb
     * where M_IjAb = 2 * gbar_IjAb - gbar_JiAb
     */
    dpdbuf4 L, M, temp;

    dct_timer_on("DCTSolver::2 * g_IJAB - g_JIAB");
    // M_IjAb = g_IjAb - g_JiAb
    global_dpd_->buf4_init(&M, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 1,
                           "MO Ints <OO|VV>");
    global_dpd_->buf4_copy(&M, PSIF_LIBTRANS_DPD, "MO Ints Temp <OO|VV>");
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_init(&M, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "MO Ints Temp <OO|VV>");
    global_dpd_->buf4_init(&temp, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "MO Ints <OO|VV>");
    // M_IjAb += g_IjAb
    dpd_buf4_add(&M, &temp, 1.0);
    global_dpd_->buf4_close(&temp);
    dct_timer_off("DCTSolver::2 * g_IJAB - g_JIAB");

    dct_timer_on("DCTSolver::lambda_IjAb M_IjAb");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");

    // E_MP2 = lambda_IjAb * M_IjAb
    double eMP2 = global_dpd_->buf4_dot(&L, &M);
    global_dpd_->buf4_close(&M);
    global_dpd_->buf4_close(&L);
    dct_timer_off("DCTSolver::lambda_IjAb M_IjAb");

    new_total_energy_ = scf_energy_ + eMP2;
    outfile->Printf("\t*Total Hartree-Fock energy        = %20.15f\n", scf_energy_);
    outfile->Printf("\t Total MP2 correlation energy     = %20.15f\n", eMP2);
    outfile->Printf("\t*Total MP2 energy                 = %20.15f\n", new_total_energy_);

    variables_["MP2 TOTAL ENERGY"] = new_total_energy_;
    variables_["MP2 CORRELATION ENERGY"] = eMP2;

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    dct_timer_off("DCTSolver::initialize_amplitudes()");
}
}  // namespace dct
}  // namespace psi
