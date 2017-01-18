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
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libpsio/psio.hpp"
#include "defines.h"



namespace psi{ namespace dcft{

/**
 * Computes the DCFT density matrix and energy
 */
double
DCFTSolver::compute_energy()
{
    double total_energy = 0.0;

    if (options_.get_str("SCF_TYPE") == "DF" || options_.get_str("SCF_TYPE") == "CD" || options_.get_str("SCF_TYPE") == "DIRECT"){
        if (!options_["DCFT_TYPE"].has_changed())
            options_.set_global_str("DCFT_TYPE", "DF");
        else if (options_.get_str("DCFT_TYPE") == "CONV")
            throw PSIEXCEPTION("Please set SCF_TYPE to PK or OUT_OF_CORE in order to use DCFT_TYPE=CONV.");
    }

    if (options_.get_str("DCFT_TYPE") == "DF"){
        if (!options_["AO_BASIS"].has_changed())
            options_.set_str("DCFT", "AO_BASIS", "NONE");
        else if (options_.get_str("AO_BASIS") == "DISK"){
            outfile->Printf("\n\n\t**** Warning: AO_BASIS=DISK not implemented in density-fitted DCFT. Switch to AO_BASIS=NONE ****\n");
            options_.set_str("DCFT", "AO_BASIS", "NONE");
        }
    }


    if(options_.get_str("REFERENCE") == "RHF")
        total_energy = compute_energy_RHF();
    else if (options_.get_str("REFERENCE") == "UHF")
        total_energy = compute_energy_UHF();
    else if (options_.get_str("REFERENCE") == "ROHF"){
        outfile->Printf("\n\n\t**** Warning: ROHF reference, then unrestricted DCFT ****\n");
        total_energy = compute_energy_UHF();
    }
    else
        throw PSIEXCEPTION("Unknown DCFT reference.");


    // Compute the analytic gradients, if requested
    if(options_.get_str("DERTYPE") == "FIRST") {
        // Shut down the timers
        tstop();
        // Start the timers
        tstart();
        // Solve the response equations, compute relaxed OPDM and TPDM and dump them to disk
        compute_gradient();
        if (options_.get_str("REFERENCE") != "RHF") {
            // Compute TPDM trace
            compute_TPDM_trace();
        }

    }

    // Free up memory and close files
    finalize();

    return(total_energy);
}

}} // Namespaces
