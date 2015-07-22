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

#include "defines.h"

using namespace boost;

namespace psi{ namespace dcft{

/**
  * Compute DCFT energy using unrestricted HF reference
  */
double DCFTSolver::compute_energy_UHF()
{
    orbitalsDone_    = false;
    cumulantDone_ = false;
    densityConverged_ = false;
    energyConverged_ = false;
    // Perform SCF guess for the orbitals
    scf_guess();
    // Perform MP2 guess for the cumulant
    mp2_guess();

    // Print out information about the job
    outfile->Printf( "\n\tDCFT Functional:    \t\t %s", options_.get_str("DCFT_FUNCTIONAL").c_str());
    outfile->Printf( "\n\tAlgorithm:          \t\t %s", options_.get_str("ALGORITHM").c_str());
    outfile->Printf( "\n\tAO-Basis Integrals: \t\t %s", options_.get_str("AO_BASIS").c_str());
    if (options_.get_str("ALGORITHM") == "QC") {
        outfile->Printf( "\n\tQC type:            \t\t %s", options_.get_str("QC_TYPE").c_str());
        outfile->Printf( "\n\tQC coupling:        \t\t %s", options_.get_bool("QC_COUPLING") ? "TRUE" : "FALSE");
    }
    if (energy_level_shift_ > 1E-6) {
        outfile->Printf( "\n\tUsing level shift of %5.3f a.u.            ", energy_level_shift_);
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

    // Orbital-optimized stuff
    if (options_.get_str("ALGORITHM") == "TWOSTEP" && orbital_optimized_)
        throw PSIEXCEPTION("Two-step algorithm cannot be run for the orbital-optimized DCFT methods");

    // Choose a paricular algorithm and solve the equations
    if(options_.get_str("ALGORITHM") == "TWOSTEP") {
        run_twostep_dcft();
    }
    else if (options_.get_str("ALGORITHM") == "SIMULTANEOUS") {
        if (!orbital_optimized_) {
            run_simult_dcft();
        }
        else {
            run_simult_dcft_oo();
        }
    }
    else if (options_.get_str("ALGORITHM") == "QC") {
        run_qc_dcft();
    }
    else {
        throw PSIEXCEPTION("Unknown DCFT algoritm");
    }

    // If not converged -> Break
    if(!orbitalsDone_ || !cumulantDone_ || !densityConverged_)
        throw ConvergenceError<int>("DCFT", maxiter_, cumulant_threshold_, cumulant_convergence_, __FILE__, __LINE__);

    outfile->Printf("\n\t*%6s SCF Energy                                 = %20.15f\n", options_.get_str("DCFT_FUNCTIONAL").c_str(), scf_energy_);
    outfile->Printf("\t*%6s Lambda Energy                              = %20.15f\n", options_.get_str("DCFT_FUNCTIONAL").c_str(), lambda_energy_);
    outfile->Printf("\t*%6s Total Energy                               = %20.15f\n", options_.get_str("DCFT_FUNCTIONAL").c_str(), new_total_energy_);


    Process::environment.globals["DCFT SCF ENERGY"]    = scf_energy_;
    Process::environment.globals["DCFT LAMBDA ENERGY"] = lambda_energy_;
    Process::environment.globals["DCFT TOTAL ENERGY"]  = new_total_energy_;

    // Compute three-particle contribution to the DCFT energy
    if (options_.get_str("THREE_PARTICLE") == "PERTURBATIVE") {
        // Check options
        if (options_.get_str("DERTYPE") == "FIRST")
            throw FeatureNotImplemented("DCFT three-particle energy correction", "Analytic gradients", __FILE__, __LINE__);
        // Compute the three-particle energy
        double three_particle_energy = compute_three_particle_energy();
        outfile->Printf("\t*DCFT Three-particle Energy                        = %20.15f\n", three_particle_energy);
        outfile->Printf("\t*DCFT Total Energy                                 = %20.15f\n", new_total_energy_ + three_particle_energy);
        // Set global variables
        Process::environment.globals["DCFT THREE-PARTICLE ENERGY"] = three_particle_energy;
        Process::environment.globals["CURRENT ENERGY"]             = new_total_energy_ + three_particle_energy;
    }
    else {
        Process::environment.globals["CURRENT ENERGY"]             = new_total_energy_;
    }

    if(!options_.get_bool("MO_RELAX")){
        outfile->Printf( "Warning!  The orbitals were not relaxed\n");
    }

    // Print natural occupations
    print_opdm();

    if (orbital_optimized_) {
        // Compute one-electron properties
        compute_oe_properties();
        // Write to MOLDEN file if requested
        if (options_.get_bool("MOLDEN_WRITE")) write_molden_file();
    }

    if(options_.get_bool("TPDM")) dump_density();
//    check_n_representability();

    if (options_.get_str("DCFT_FUNCTIONAL") == "CEPA0") {
        compute_unrelaxed_density_OOOO();
        compute_unrelaxed_density_OVOV();
        compute_unrelaxed_density_VVVV();
        compute_TPDM_trace();
    }

    return(new_total_energy_);
}

}} // Namespace

