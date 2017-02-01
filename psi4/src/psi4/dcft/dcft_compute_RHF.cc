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
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"
#include "defines.h"



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

    // If DCFT computation type is density fitting, build b(Q|mn) in AO and SO basis
    if (options_.get_str("DCFT_TYPE") == "DF") {
        df_build_b_ao();
        transform_b_ao2so();
        transform_b_ao2so_scf();
    }

    // Perform MP2 guess for the cumulant
    mp2_guess_RHF();

    // Print out information about the job
    outfile->Printf( "\n\tDCFT Functional:    \t\t %s", options_.get_str("DCFT_FUNCTIONAL").c_str());
    outfile->Printf( "\n\tDCFT Type:          \t\t %s", options_.get_str("DCFT_TYPE").c_str());
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
    if (options_.get_str("ALGORITHM") == "QC" || options_.get_str("ALGORITHM") == "TWOSTEP")
        throw FeatureNotImplemented("RHF-reference DCFT", "ALGORITHM = QC/TWOSTEP", __FILE__, __LINE__);
    if (options_.get_str("DCFT_FUNCTIONAL") == "ODC-13")
        throw FeatureNotImplemented("RHF-reference DCFT", "DCFT_FUNCTIONAL = ODC-13", __FILE__, __LINE__);
    if (options_.get_str("THREE_PARTICLE") == "PERTURBATIVE")
        throw FeatureNotImplemented("RHF-reference DCFT", "Three-particle energy correction", __FILE__, __LINE__);

    // Orbital-optimized stuff
    if (options_.get_str("ALGORITHM") == "TWOSTEP" && orbital_optimized_)
        throw PSIEXCEPTION("Two-step algorithm cannot be run for orbital-optimized DCFT methods");

    // Choose a paricular algorithm and solve the equations
    if (options_.get_str("ALGORITHM") == "SIMULTANEOUS") {
        if (!orbital_optimized_) {
            run_simult_dcft_RHF();
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

    std::string prefix = options_.get_str("DCFT_TYPE") == "DF"? "DF-" : " ";

    outfile->Printf("\n\t*%3s%5s SCF Energy                                 = %23.15f\n", prefix.c_str(), options_.get_str("DCFT_FUNCTIONAL").c_str(), scf_energy_);
    outfile->Printf(  "\t*%3s%5s Lambda Energy                              = %23.15f\n", prefix.c_str(), options_.get_str("DCFT_FUNCTIONAL").c_str(), lambda_energy_);
    outfile->Printf(  "\t*%3s%5s Total Energy                               = %23.15f\n", prefix.c_str(), options_.get_str("DCFT_FUNCTIONAL").c_str(), new_total_energy_);

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

void DCFTSolver::run_simult_dcft_RHF()
{
    // This is the simultaneous orbital/lambda update algorithm
    int cycle = 0;
    outfile->Printf( "\n\n\t*=================================================================================*\n"
                         "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
                         "\t*---------------------------------------------------------------------------------*\n");
    SharedMatrix tmp = SharedMatrix(new Matrix("temp", nirrep_, nsopi_, nsopi_));
    // Set up the DIIS manager
    DIISManager diisManager(maxdiis_, "DCFT DIIS vectors");

    // DIIS on orbitals (AA and BB) and cumulants (AA, AB, BB)
    dpdbuf4 Laa, Lab, Lbb;
    global_dpd_->buf4_init(&Laa, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
    global_dpd_->buf4_init(&Lbb, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <oo|vv>");
    diisManager.set_error_vector_size(5, DIISEntry::Matrix, scf_error_a_.get(),
                                         DIISEntry::Matrix, scf_error_b_.get(),
                                         DIISEntry::DPDBuf4, &Laa,
                                         DIISEntry::DPDBuf4, &Lab,
                                         DIISEntry::DPDBuf4, &Lbb);
    diisManager.set_vector_size(5, DIISEntry::Matrix, Fa_.get(),
                                   DIISEntry::Matrix, Fb_.get(),
                                   DIISEntry::DPDBuf4, &Laa,
                                   DIISEntry::DPDBuf4, &Lab,
                                   DIISEntry::DPDBuf4, &Lbb);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);

    while((!orbitalsDone_ || !cumulantDone_ || !densityConverged_ || !energyConverged_)
            && cycle++ < maxiter_){
        std::string diisString;
        // Save the old energy
        old_total_energy_ = new_total_energy_;
        // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
        build_tau_RHF();
        if (exact_tau_) {
            refine_tau_RHF();
        }
        transform_tau_RHF();

        if(options_.get_str("DCFT_TYPE") == "DF" && options_.get_str("AO_BASIS") == "NONE"){
            build_DF_tensors_RHF();

            SharedMatrix mo_h = SharedMatrix(new Matrix("MO-based H", nirrep_, nmopi_, nmopi_));
            mo_h->copy(so_h_);
            mo_h->transform(Ca_);

            moFa_->copy(mo_h);
            moFa_->add(mo_gbarGamma_A_);
            // Back-transform the Fock matrix to the SO basis: F_so = (Ct)^-1 F_mo C^-1 = (C^-1)t F_mo C^-1
            SharedMatrix Ca_inverse = SharedMatrix(new Matrix("Ca_ inverse", nirrep_, nmopi_, nsopi_));
            Ca_inverse->copy(Ca_);
            Ca_inverse->general_invert();
            Fa_->copy(moFa_);
            Fa_->transform(Ca_inverse);
            Fb_->copy(Fa_);
        }
        else{
            // Copy core hamiltonian into the Fock matrix array: F = H
            Fa_->copy(so_h_);

            // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
            process_so_ints_RHF();

            // Add non-idempotent density contribution (Tau) to the Fock matrix: F += Gbar * Tau
            Fa_->add(g_tau_a_);
            // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
            moFa_->copy(Fa_);
            // Transform the Fock matrix to the MO basis
            moFa_->transform(Ca_);
            Fb_->copy(Fa_);
        }
        // Compute new SCF energy
        compute_scf_energy_RHF();
        // Add SCF energy contribution to the total DCFT energy
        new_total_energy_ = scf_energy_;
        // Check SCF convergence
        orbitals_convergence_ = compute_scf_error_vector_RHF();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
        build_cumulant_intermediates_RHF();
        // Compute the residuals for density cumulant equations
        cumulant_convergence_ = compute_cumulant_residual_RHF();
        if (fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCFT density cumulant equations diverged");
        // Check convergence for density cumulant iterations
        cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
        // Update density cumulant tensor
        update_cumulant_jacobi_RHF();
        // Compute new DCFT energy (lambda contribution)
        compute_dcft_energy_RHF();
        // Add lambda energy to the DCFT total energy
        new_total_energy_ += lambda_energy_;
        // Check convergence of the total DCFT energy
        energyConverged_ = fabs(old_total_energy_ - new_total_energy_) < energy_threshold_;
        if(orbitals_convergence_ < diis_start_thresh_ && cumulant_convergence_ < diis_start_thresh_){
            //Store the DIIS vectors
            dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
            // Compute R_OOVV and R_oovv from R_OoVv, used as DIIS error vectors
            compute_R_AA_and_BB();

            global_dpd_->buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                          ID("[O,O]"), ID("[V,V]"), 0, "R <OO|VV>");
            global_dpd_->buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                          ID("[O,O]"), ID("[V,V]"), 0, "R SF <OO|VV>"); // R <Oo|Vv>
            global_dpd_->buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                          ID("[O,O]"), ID("[V,V]"), 0, "R <oo|vv>");
            global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                          ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
            global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                          ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>"); // Lambda <Oo|Vv>
            global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                          ID("[O,O]"), ID("[V,V]"), 0, "Lambda <oo|vv>");
            if(diisManager.add_entry(10, scf_error_a_.get(), scf_error_b_.get(), &Raa, &Rab, &Rbb,
                                       Fa_.get(), Fb_.get(), &Laa, &Lab, &Lbb)){
                diisString += "S";
            }
            if(diisManager.subspace_size() > mindiisvecs_){
                diisString += "/E";
                diisManager.extrapolate(5, Fa_.get(), Fb_.get(), &Laa, &Lab, &Lbb);
            }
            global_dpd_->buf4_close(&Raa);
            global_dpd_->buf4_close(&Rab);
            global_dpd_->buf4_close(&Rbb);
            global_dpd_->buf4_close(&Laa);
            global_dpd_->buf4_close(&Lab);
            global_dpd_->buf4_close(&Lbb);
        }
        // Transform the Fock matrix to the symmetrically orhogonalized basis set and digonalize it
        // Obtain new orbitals
        Fa_->transform(s_half_inv_);
        Fa_->diagonalize(tmp, epsilon_a_);
        old_ca_->copy(Ca_);
        Ca_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        old_cb_->copy(old_ca_);
        Cb_->copy(Ca_);
        // Make sure that the orbital phase is retained
        if(!correct_mo_phases(false)){
            outfile->Printf("\t\tThere was a problem correcting the MO phases.\n"
                            "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
        }
        // Transform two-electron integrals to the MO basis using new orbitals, build denominators
        transform_integrals_RHF();
        // Update SCF density (Kappa) and check its RMS
        densityConverged_ = update_scf_density_RHF() < orbitals_threshold_;
        // If we've performed enough lambda updates since the last orbitals
        // update, reset the counter so another SCF update is performed
        outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                cycle, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                new_total_energy_, diisString.c_str());
    }
    outfile->Printf( "\t*=================================================================================*\n");
}

}} // Namespace
