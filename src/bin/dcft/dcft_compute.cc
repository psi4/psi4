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

#include "dcft.h"
#include <cmath>
#include <libdpd/dpd.h>
#include <libtrans/integraltransform.h>
#include <libdiis/diismanager.h>
#include <libpsio/psio.hpp>
#include "defines.h"

using namespace boost;

namespace psi{ namespace dcft{

/**
 * Computes the DCFT density matrix and energy
 */
double
DCFTSolver::compute_energy()
{
    double total_energy = 0.0;

    if(options_.get_str("REFERENCE") == "RHF")
        total_energy = compute_energy_RHF();
    else if (options_.get_str("REFERENCE") == "UHF")
        total_energy = compute_energy_UHF();
    else if (options_.get_str("REFERENCE") == "ROHF")
        throw FeatureNotImplemented("DCFT", "ROHF reference", __FILE__, __LINE__);
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
        if (options_.get_str("REFERENCE") == "UHF") {
            // Compute TPDM trace
            compute_TPDM_trace();
        }

    }

    // Free up memory and close files
    finalize();

    return(total_energy);
}

}} // Namespaces

