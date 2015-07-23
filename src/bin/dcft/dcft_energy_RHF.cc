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
#include <libpsio/psio.h>
#include <psifiles.h>
#include "defines.h"

using namespace boost;

namespace psi{ namespace dcft{

/**
 * Compute DCFT energy using restricted HF reference
 */
double DCFTSolver::compute_energy_RHF()
{
    orbitalsDone_    = false;
    cumulantDone_ = false;
    densityConverged_ = false;
    energyConverged_ = false;

    // Perform SCF guess for the orbitals
    scf_guess_RHF();

    // Perform MP2 guess for the cumulant
    mp2_guess_RHF();

    // Print out information about the job
    outfile->Printf( "\n\tDCFT Functional:    \t\t %s", options_.get_str("DCFT_FUNCTIONAL").c_str());
    outfile->Printf( "\n\tAlgorithm:          \t\t %s", options_.get_str("ALGORITHM").c_str());
    outfile->Printf( "\n\tAO-Basis Integrals: \t\t %s", options_.get_str("AO_BASIS").c_str());
    if (options_.get_str("ALGORITHM") == "QC") {
        outfile->Printf( "\n\tQC type:            \t\t %s", options_.get_str("QC_TYPE").c_str());
        outfile->Printf( "\n\tQC coupling:        \t\t %s", options_.get_bool("QC_COUPLING") ? "TRUE" : "FALSE");
    }

    // Things that are not implemented yet...
    if (options_.get_str("DERTYPE") == "FIRST" && (options_.get_str("DCFT_FUNCTIONAL") == "DC-12"))
        throw FeatureNotImplemented("DC-12 functional", "Analytic gradients", __FILE__, __LINE__);
    if (options_.get_str("AO_BASIS") == "DISK" && options_.get_str("DCFT_FUNCTIONAL") == "CEPA0")
        throw FeatureNotImplemented("CEPA0", "AO_BASIS = DISK", __FILE__, __LINE__);
    if (options_.get_str("AO_BASIS") == "DISK" && options_.get_str("ALGORITHM") == "QC" && options_.get_str("QC_TYPE") == "SIMULTANEOUS")
        throw FeatureNotImplemented("Simultaneous QC", "AO_BASIS = DISK", __FILE__, __LINE__);
    if (!(options_.get_str("ALGORITHM") == "TWOSTEP") && options_.get_str("DCFT_FUNCTIONAL") == "CEPA0")
        throw FeatureNotImplemented("CEPA0", "Requested DCFT algorithm", __FILE__, __LINE__);
    if (options_.get_str("ALGORITHM") == "QC" || options_.get_str("ALGORITHM") == "TWOSTEP")
        throw FeatureNotImplemented("Spin-adapted RHF-reference DCFT", "ALGORITHM = QC/TWOSTEP", __FILE__, __LINE__);

    // Orbital-optimized stuff
    if (options_.get_str("ALGORITHM") == "TWOSTEP" && orbital_optimized_)
        throw PSIEXCEPTION("Two-step algorithm cannot be run for orbital-optimized DCFT methods");

    // Choose a paricular algorithm and solve the equations
    if (options_.get_str("ALGORITHM") == "SIMULTANEOUS") {
        if (!orbital_optimized_) {
            throw PSIEXCEPTION("RHF reference only available for orbital-optimized DCFT methods.");
        }
        else {
            run_simult_dcft_oo_RHF();
        }
    }
    else {
        throw PSIEXCEPTION("Unknown DCFT algoritm");
    }

    // If not converged -> Break
    if(!orbitalsDone_ || !cumulantDone_ || !densityConverged_)
        throw ConvergenceError<int>("DCFT", maxiter_, cumulant_threshold_, cumulant_convergence_, __FILE__, __LINE__);

    outfile->Printf( "\n\t*DCFT SCF Energy                                 = %20.15f\n", scf_energy_);
    outfile->Printf(   "\t*DCFT Lambda Energy                              = %20.15f\n", lambda_energy_);
    outfile->Printf(   "\t*DCFT Total Energy                               = %20.15f\n", new_total_energy_);

    Process::environment.globals["CURRENT ENERGY"] = new_total_energy_;
    Process::environment.globals["DCFT TOTAL ENERGY"] = new_total_energy_;
    Process::environment.globals["DCFT SCF ENERGY"] = scf_energy_;
    Process::environment.globals["DCFT LAMBDA ENERGY"] = lambda_energy_;

    if(!options_.get_bool("MO_RELAX")){
        outfile->Printf( "Warning!  The orbitals were not relaxed\n");
    }

    print_opdm_RHF();

    return new_total_energy_;
}

/**
 * Uses the intermediates to compute the energy
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

    dcft_timer_on("Timing::G_IjAb + g_IjAb");
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
    dcft_timer_off("Timing::G_IjAb + g_IjAb");

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

}} // Namespace
