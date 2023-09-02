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

#include "dct.h"
#include "psi4/psifiles.h"

#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

#include <cmath>

namespace psi {
namespace dct {

/**
 * Compute DCT energy using restricted HF reference
 */
double DCTSolver::compute_energy_RHF() {
    orbitalsDone_ = false;
    cumulantDone_ = false;
    energyConverged_ = false;
    initialize_orbitals_from_reference_R();

    // If DCT computation type is density fitting, build b(Q|mn)
    if (options_.get_str("DCT_TYPE") == "DF") {
        initialize_df();
        build_df_b();
    }

    initialize_integraltransform();
    // Perform MP2 guess for the cumulant
    initialize_amplitudes_RHF();

    // Print out information about the job
    outfile->Printf("\n\tDCT Functional:    \t\t %s", options_.get_str("DCT_FUNCTIONAL").c_str());
    outfile->Printf("\n\tDCT Type:          \t\t %s", options_.get_str("DCT_TYPE").c_str());
    outfile->Printf("\n\tAlgorithm:          \t\t %s", options_.get_str("ALGORITHM").c_str());
    outfile->Printf("\n\tAO-Basis Integrals: \t\t %s", options_.get_str("AO_BASIS").c_str());
    if (options_.get_str("ALGORITHM") == "QC") {
        outfile->Printf("\n\tQC type:            \t\t %s", options_.get_str("QC_TYPE").c_str());
        outfile->Printf("\n\tQC coupling:        \t\t %s", options_.get_bool("QC_COUPLING") ? "TRUE" : "FALSE");
    }
    if (energy_level_shift_ > 1E-6) {
        outfile->Printf("\n\tUsing level shift of %5.3f a.u.            ", energy_level_shift_);
    }

    // Orbital-optimized stuff
    if (options_.get_str("ALGORITHM") == "TWOSTEP" && orbital_optimized_)
        throw PSIEXCEPTION("Two-step algorithm cannot be run for orbital-optimized DCT methods");

    // Choose a paricular algorithm and solve the equations
    if (options_.get_str("ALGORITHM") == "SIMULTANEOUS") {
        if (!orbital_optimized_) {
            run_simult_dct_RHF();
        } else {
            run_simult_dct_oo_RHF();
        }
    } else {
        throw PSIEXCEPTION("Unknown DCT algorithm");
    }

    // If not converged -> Break
    if (!orbitalsDone_ || !cumulantDone_)
        throw ConvergenceError<int>("DCT", maxiter_, cumulant_threshold_, cumulant_convergence_, __FILE__, __LINE__);

    std::string prefix = options_.get_str("DCT_TYPE") == "DF" ? "DF-" : " ";

    outfile->Printf("\n\t*%3s%5s SCF Energy                                 = %23.15f\n", prefix.c_str(),
                    options_.get_str("DCT_FUNCTIONAL").c_str(), scf_energy_);
    outfile->Printf("\t*%3s%5s Lambda Energy                              = %23.15f\n", prefix.c_str(),
                    options_.get_str("DCT_FUNCTIONAL").c_str(), lambda_energy_);
    outfile->Printf("\t*%3s%5s Total Energy                               = %23.15f\n", prefix.c_str(),
                    options_.get_str("DCT_FUNCTIONAL").c_str(), new_total_energy_);

    set_scalar_variable("CURRENT ENERGY", new_total_energy_);
    set_scalar_variable("DCT TOTAL ENERGY", new_total_energy_);
    set_scalar_variable("DCT SCF ENERGY", scf_energy_);
    set_scalar_variable("DCT LAMBDA ENERGY", lambda_energy_);

    print_opdm_RHF();

    if (orbital_optimized_) {
        construct_oo_density_RHF();
    }

    return new_total_energy_;
}

void DCTSolver::run_simult_dct_RHF() {
    // This is the simultaneous orbital/lambda update algorithm
    int cycle = 0;
    outfile->Printf(
        "\n\n\t*=================================================================================*\n"
        "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
        "\t*---------------------------------------------------------------------------------*\n");
    auto tmp = std::make_shared<Matrix>("temp", nirrep_, nsopi_, nsopi_);
    // Set up the DIIS manager
    DIISManager diisManager(maxdiis_, "DCT DIIS vectors");

    // DIIS on orbitals (AA and BB) and cumulants (AA, AB, BB)
    dpdbuf4 Laa, Lab, Lbb;
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude SF <OO|VV>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                           "Amplitude <oo|vv>");
    diisManager.set_error_vector_size(scf_error_a_.get(), scf_error_b_.get(), &Laa, &Lab, &Lbb);
    diisManager.set_vector_size(Fa_.get(), Fb_.get(), &Laa, &Lab, &Lbb);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);

    while ((!orbitalsDone_ || !cumulantDone_ || !energyConverged_) && cycle++ < maxiter_) {
        std::string diisString;
        // Save the old energy
        old_total_energy_ = new_total_energy_;
        // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
        compute_SO_tau_R();

        if (options_.get_str("DCT_TYPE") == "DF" && options_.get_str("AO_BASIS") == "NONE") {
            build_DF_tensors_RHF();

            auto mo_h = so_h_.clone();
            mo_h->transform(Ca_);

            moFa_->copy(mo_h);
            moFa_->add(mo_gbarGamma_A_);
            // Back-transform the Fock matrix to the SO basis: F_so = (Ct)^-1 F_mo C^-1 = (C^-1)t F_mo C^-1
            auto Ca_inverse = Ca_->clone();
            Ca_inverse->general_invert();
            Fa_->copy(moFa_);
            Fa_->transform(Ca_inverse);
            Fb_->copy(Fa_);
        } else {
            // Copy core hamiltonian into the Fock matrix array: F = H
            Fa_->copy(so_h_);

            // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
            process_so_ints_RHF();

            // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
            moFa_->copy(Fa_);
            // Transform the Fock matrix to the MO basis
            moFa_->transform(Ca_);
            Fb_->copy(Fa_);
        }
        // Compute new SCF energy
        compute_scf_energy_RHF();
        // Add SCF energy contribution to the total DCT energy
        new_total_energy_ = scf_energy_;
        // Check SCF convergence
        orbitals_convergence_ = compute_scf_error_vector_RHF();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        // Build G and F intermediates needed for the density cumulant residual equations and DCT energy computation
        build_cumulant_intermediates_RHF();
        // Compute the residuals for density cumulant equations
        cumulant_convergence_ = compute_cumulant_residual_RHF();
        if (std::fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCT density cumulant equations diverged");
        // Check convergence for density cumulant iterations
        cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
        // Update density cumulant tensor
        update_cumulant_jacobi_RHF();
        // Compute new DCT energy (lambda contribution)
        compute_dct_energy_RHF();
        // Add lambda energy to the DCT total energy
        new_total_energy_ += lambda_energy_;
        // Check convergence of the total DCT energy
        energyConverged_ = std::fabs(old_total_energy_ - new_total_energy_) < energy_threshold_;
        if (orbitals_convergence_ < diis_start_thresh_ && cumulant_convergence_ < diis_start_thresh_) {
            // Store the DIIS vectors
            dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
            // Compute R_OOVV and R_oovv from R_OoVv, used as DIIS error vectors
            compute_R_AA_and_BB();

            global_dpd_->buf4_init(&Raa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "R <OO|VV>");
            global_dpd_->buf4_init(&Rab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "R SF <OO|VV>");  // R <Oo|Vv>
            global_dpd_->buf4_init(&Rbb, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "R <oo|vv>");
            global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "Amplitude <OO|VV>");
            global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "Amplitude SF <OO|VV>");  // Amplitude <Oo|Vv>
            global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "Amplitude <oo|vv>");
            if (diisManager.add_entry(scf_error_a_.get(), scf_error_b_.get(), &Raa, &Rab, &Rbb, Fa_.get(),
                                      Fb_.get(), &Laa, &Lab, &Lbb)) {
                diisString += "S";
            }
            if (diisManager.subspace_size() > mindiisvecs_) {
                diisString += "/E";
                diisManager.extrapolate(Fa_.get(), Fb_.get(), &Laa, &Lab, &Lbb);
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
        if (!correct_mo_phases(false)) {
            outfile->Printf(
                "\t\tThere was a problem correcting the MO phases.\n"
                "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
        }
        // Transform two-electron integrals to the MO basis using new orbitals, build denominators
        transform_integrals_RHF();
        // Update SCF density (Kappa)
        update_scf_density_RHF();
        // If we've performed enough lambda updates since the last orbitals
        // update, reset the counter so another SCF update is performed
        outfile->Printf("\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n", cycle, orbitals_convergence_,
                        cumulant_convergence_, new_total_energy_ - old_total_energy_, new_total_energy_,
                        diisString.c_str());
    }
    outfile->Printf("\t*=================================================================================*\n");
}

}  // namespace dct
}  // namespace psi
