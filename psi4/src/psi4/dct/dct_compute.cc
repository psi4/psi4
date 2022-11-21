/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "psi4/libdpd/dpd.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"

#include <cmath>

namespace psi {
namespace dct {

/**
 * Computes the DCT density matrix and energy
 */
double DCTSolver::compute_energy() {
    double total_energy = 0.0;

    if (options_.get_str("DCT_GUESS") == "DCT") {
        // We're reusing DCT files from a previous computation. Let's check that one exists.
        if (!((same_a_b_orbs_ == true && psio_tocentry_exists(PSIF_DCT_DPD, "Amplitude SF <OO|VV>")) ||
              (same_a_b_orbs_ == false && psio_tocentry_exists(PSIF_DCT_DPD, "Amplitude <Oo|Vv>")))) {
            throw PSIEXCEPTION("Could not find a previous DCT computation as the DCT_GUESS=DCT guess.");
        }
    } else {
        // Obliviate all old DCT files.
        psio_->close(PSIF_DCT_DPD, 0);
        psio_->open(PSIF_DCT_DPD, PSIO_OPEN_OLD);
    }

    validate_energy();

    if (same_a_b_orbs_ == true)
        total_energy = compute_energy_RHF();
    else {
        if (reference_wavefunction_->name() == "ROHF")
            outfile->Printf("\n\n\t**** Warning: ROHF reference, then unrestricted DCT ****\n");
        total_energy = compute_energy_UHF();
    }

    // If we have a variational method, the density is free.
    if (options_.get_bool("OPDM") == true || (orbital_optimized_ && options_.get_str("THREE_PARTICLE") == "NONE")) {
        tstop();
        tstart();
        compute_relaxed_opdm();
    }

    // Compute the analytic gradients, if requested
    if (options_.get_str("DERTYPE") == "FIRST") {
        // Shut down the timers
        tstop();
        // Start the timers
        tstart();
        // Solve the response equations, compute relaxed OPDM and TPDM and dump them to disk
        compute_gradient();
        if (options_.get_str("REFERENCE") != "RHF") {
            // Compute TPDM trace
            compute_TPDM_trace(options_.get_str("DCT_TYPE") == "DF");
        }
    }

    // Enforce that ODPM is hermitian. This may fail when the non-OO DCT is used.
    if (Da_) Da_->hermitivitize();
    if (Db_) Db_->hermitivitize();

    // Free up memory and close files
    finalize();

    return (total_energy);
}

void DCTSolver::compute_relaxed_opdm() {
    validate_opdm();
    if (orbital_optimized_ && options_.get_str("THREE_PARTICLE") == "NONE") {
        // Our energy functional variationally minimizes the orbitals and amplitudes.
        // No relaxation terms to compute.
        Da_ = std::make_shared<Matrix>(std::move(construct_oo_density(aocc_tau_, avir_tau_, *kappa_mo_a_, *Ca_)));
        if (same_a_b_orbs_) {
            Db_ = Da_;
        } else {
            Db_ = std::make_shared<Matrix>(std::move(construct_oo_density(bocc_tau_, bvir_tau_, *kappa_mo_b_, *Cb_)));
        }
        return;
    }
    // We must be in the DC-06 case.
    dc06_response_init();
    dc06_response();
    dc06_compute_relaxed_density_1PDM();
}

}  // namespace dct
}  // namespace psi
