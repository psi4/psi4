/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"

#include <cmath>

namespace psi {
namespace dct {

/**
 * Compute DCT energy using unrestricted HF reference
 */
double DCTSolver::compute_energy_UHF() {
    orbitalsDone_ = false;
    cumulantDone_ = false;
    densityConverged_ = false;
    energyConverged_ = false;
    // Perform SCF guess for the orbitals
    scf_guess();

    // If DCT computation type is density fitting, build b(Q|mn)
    if (options_.get_str("DCT_TYPE") == "DF") {
        df_build_b();
    }

    // Perform MP2 guess for the cumulant
    mp2_guess();

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

    // Things that are not implemented yet...
    if (options_.get_str("DERTYPE") == "FIRST" && (options_.get_str("DCT_FUNCTIONAL") == "DC-12"))
        throw FeatureNotImplemented("DC-12 functional", "Analytic gradients", __FILE__, __LINE__);
    if (!(options_.get_str("DCT_FUNCTIONAL") == "ODC-06" || options_.get_str("DCT_FUNCTIONAL") == "ODC-12" ||
          options_.get_str("DCT_FUNCTIONAL") == "DC-06" || options_.get_str("DCT_FUNCTIONAL") == "DC-12") &&
        options_.get_str("DCT_TYPE") == "DF")
        throw FeatureNotImplemented("ODC-13", "Density Fitting", __FILE__, __LINE__);
    if (options_.get_str("THREE_PARTICLE") == "PERTURBATIVE" && options_.get_str("DCT_TYPE") == "DF")
        throw FeatureNotImplemented("Three-particle energy correction", "Density Fitting", __FILE__, __LINE__);
    if (options_.get_str("ALGORITHM") == "QC") {
        if (options_.get_str("AO_BASIS") == "DISK" && options_.get_str("QC_TYPE") == "SIMULTANEOUS")
            throw FeatureNotImplemented("Simultaneous QC", "AO_BASIS = DISK", __FILE__, __LINE__);
        if (options_.get_str("DCT_TYPE") == "DF")
            throw FeatureNotImplemented("QC algorithm", "Density Fitting", __FILE__, __LINE__);
        if (options_.get_str("DCT_FUNCTIONAL") != "DC-06")
            outfile->Printf(
                "\n\n\t**** Warning: Using DC-06 hessian, as others not implemented. Quadratic convergence is not "
                "guaranteed. ****\n");
    }

    // Orbital-optimized stuff
    if (options_.get_str("ALGORITHM") == "TWOSTEP" && orbital_optimized_)
        throw PSIEXCEPTION("Two-step algorithm cannot be run for the orbital-optimized DCT methods");

    // Choose a paricular algorithm and solve the equations
    if (options_.get_str("ALGORITHM") == "TWOSTEP") {
        run_twostep_dct();
    } else if (options_.get_str("ALGORITHM") == "SIMULTANEOUS") {
        if (!orbital_optimized_) {
            run_simult_dct();
        } else {
            run_simult_dct_oo();
        }
    } else if (options_.get_str("ALGORITHM") == "QC") {
        run_qc_dct();
    } else {
        throw PSIEXCEPTION("Unknown DCT algorithm");
    }

    // If not converged -> Break
    if (!orbitalsDone_ || !cumulantDone_ || !densityConverged_)
        throw ConvergenceError<int>("DCT", maxiter_, cumulant_threshold_, cumulant_convergence_, __FILE__, __LINE__);

    std::string prefix = options_.get_str("DCT_TYPE") == "DF" ? "DF-" : " ";

    outfile->Printf("\n\t*%3s%5s SCF Energy                                 = %23.15f\n", prefix.c_str(),
                    options_.get_str("DCT_FUNCTIONAL").c_str(), scf_energy_);
    outfile->Printf("\t*%3s%5s Lambda Energy                              = %23.15f\n", prefix.c_str(),
                    options_.get_str("DCT_FUNCTIONAL").c_str(), lambda_energy_);
    outfile->Printf("\t*%3s%5s Total Energy                               = %23.15f\n", prefix.c_str(),
                    options_.get_str("DCT_FUNCTIONAL").c_str(), new_total_energy_);

    Process::environment.globals["DCT SCF ENERGY"] = scf_energy_;
    Process::environment.globals["DCT LAMBDA ENERGY"] = lambda_energy_;
    Process::environment.globals["DCT TOTAL ENERGY"] = new_total_energy_;

    // Compute three-particle contribution to the DCT energy
    if (options_.get_str("THREE_PARTICLE") == "PERTURBATIVE") {
        // Check options
        if (options_.get_str("DERTYPE") == "FIRST")
            throw FeatureNotImplemented("DCT three-particle energy correction", "Analytic gradients", __FILE__,
                                        __LINE__);
        // Compute the three-particle energy
        double three_particle_energy = compute_three_particle_energy();
        outfile->Printf("\t*DCT Three-particle Energy                        = %20.15f\n", three_particle_energy);
        outfile->Printf("\t*DCT Total Energy                                 = %20.15f\n",
                        new_total_energy_ + three_particle_energy);
        // Set global variables
        Process::environment.globals["DCT THREE-PARTICLE ENERGY"] = three_particle_energy;
        Process::environment.globals["CURRENT ENERGY"] = new_total_energy_ + three_particle_energy;
    } else {
        Process::environment.globals["CURRENT ENERGY"] = new_total_energy_;
    }

    // Print natural occupations
    print_opdm();

    if (orbital_optimized_) {
        construct_oo_density_UHF();
        compute_oe_properties();
        if (options_.get_bool("MOLDEN_WRITE")) write_molden_file();
    }

    if (options_.get_bool("TPDM")) dump_density();
    //    check_n_representability();

    return (new_total_energy_);
}

void DCTSolver::run_twostep_dct() {
    // This is the two-step update - in each macro iteration, update the orbitals first, then update lambda
    // to self-consistency, until converged.  When lambda is converged and only one scf cycle is needed to reach
    // the desired cutoff, we're done

    int cycle = 0;
    outfile->Printf(
        "\n\n\t*=================================================================================*\n"
        "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
        "\t*---------------------------------------------------------------------------------*\n");

    // Set up the DIIS manager for the density cumulant and SCF iterations
    old_ca_->copy(Ca_);
    old_cb_->copy(Cb_);
    // Save F0 = H + G * Kappa for the Fock intermediate update in lambda iterations
    moF0a_->copy(Fa_);
    moF0b_->copy(Fb_);
    moF0a_->transform(Ca_);
    moF0b_->transform(Cb_);
    // Just so the correct value is printed in the first macro iteration
    orbitals_convergence_ = compute_scf_error_vector();
    // Start macro-iterations
    while ((!orbitalsDone_ || !cumulantDone_) && cycle++ < maxiter_) {
        outfile->Printf(
            "\t                          *** Macro Iteration %d ***\n"
            "\tCumulant Iterations\n",
            cycle);
        // If it's the first iteration and the user requested to relax guess orbitals, then skip the density cumulant
        // update
        if ((cycle != 1) || !options_.get_bool("RELAX_GUESS_ORBITALS")) {
            run_twostep_dct_cumulant_updates();
        } else
            outfile->Printf("\tSkipping the cumulant update to relax guess orbitals\n");
        // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
        build_tau();
        // Compute tau exactly if requested
        if (exact_tau_) {
            refine_tau();
        }
        transform_tau();
        run_twostep_dct_orbital_updates();
    }

    outfile->Printf("\t*=================================================================================*\n");
}

int DCTSolver::run_twostep_dct_cumulant_updates() {
    // Set up DIIS
    dpdbuf4 Laa, Lab, Lbb;
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Lambda <OO|VV>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Lambda <oo|vv>");
    DIISManager lambdaDiisManager(maxdiis_, "DCT DIIS Lambdas", DIISManager::LargestError, DIISManager::InCore);
    if ((nalpha_ + nbeta_) > 1) {
        lambdaDiisManager.set_error_vector_size(3, DIISEntry::DPDBuf4, &Laa, DIISEntry::DPDBuf4, &Lab,
                                                DIISEntry::DPDBuf4, &Lbb);
        lambdaDiisManager.set_vector_size(3, DIISEntry::DPDBuf4, &Laa, DIISEntry::DPDBuf4, &Lab, DIISEntry::DPDBuf4,
                                          &Lbb);
    }
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);
    cumulantDone_ = false;
    int nLambdaIterations = 0;
    // Start density cumulant (lambda) iterations
    while ((!cumulantDone_ || !energyConverged_) && nLambdaIterations++ < maxiter_) {
        std::string diisString;
        // Build new Tau from current Lambda
        if (options_.get_bool("RELAX_TAU")) {
            build_tau();
            // Compute tau exactly if requested
            if (exact_tau_) {
                refine_tau();
            }
            if (options_.get_str("AO_BASIS") == "DISK") {
                // Transform new Tau to the SO basis
                transform_tau();
                // Build SO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                build_AO_tensors();
            } else {
                // Compute GTau contribution for the Fock operator
                build_gtau();
            }
            // Update Fock operator for the F intermediate
            update_fock();
        } else {
            if (options_.get_str("AO_BASIS") == "DISK") {
                // Build SO basis tensors for the <VV||VV>, <vv||vv>, and <Vv|Vv> terms in the G intermediate
                build_AO_tensors();
            }
        }
        // Build G and F intermediates needed for the density cumulant residual equations and DCT energy computation
        build_cumulant_intermediates();
        // Compute the residuals for density cumulant equations
        cumulant_convergence_ = compute_cumulant_residual();
        // Update density cumulant tensor
        update_cumulant_jacobi();
        if (cumulant_convergence_ < diis_start_thresh_ && (nalpha_ + nbeta_) > 1) {
            // Store the DIIS vectors
            dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
            global_dpd_->buf4_init(&Raa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                                   "R <OO|VV>");
            global_dpd_->buf4_init(&Rab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                                   "R <Oo|Vv>");
            global_dpd_->buf4_init(&Rbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                                   "R <oo|vv>");
            global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                                   "Lambda <OO|VV>");
            global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                                   "Lambda <Oo|Vv>");
            global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                                   "Lambda <oo|vv>");

            if (lambdaDiisManager.add_entry(6, &Raa, &Rab, &Rbb, &Laa, &Lab, &Lbb)) {
                diisString += "S";
            }
            if (lambdaDiisManager.subspace_size() >= mindiisvecs_ && maxdiis_ > 0) {
                diisString += "/E";
                lambdaDiisManager.extrapolate(3, &Laa, &Lab, &Lbb);
            }
            global_dpd_->buf4_close(&Raa);
            global_dpd_->buf4_close(&Rab);
            global_dpd_->buf4_close(&Rbb);
            global_dpd_->buf4_close(&Laa);
            global_dpd_->buf4_close(&Lab);
            global_dpd_->buf4_close(&Lbb);
        }
        // Save old DCT energy
        old_total_energy_ = new_total_energy_;
        // Compute new DCT energy (lambda contribution)
        compute_dct_energy();
        new_total_energy_ = scf_energy_ + lambda_energy_;
        // Check convergence for density cumulant iterations
        cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
        energyConverged_ = std::fabs(new_total_energy_ - old_total_energy_) < energy_threshold_;
        if (options_.get_str("ALGORITHM") == "TWOSTEP") {
            outfile->Printf("\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n", nLambdaIterations,
                            orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                            new_total_energy_, diisString.c_str());
        }
        if (std::fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCT density cumulant equations diverged");
    }

    return nLambdaIterations;
}

void DCTSolver::run_twostep_dct_orbital_updates() {
    auto tmp = std::make_shared<Matrix>("temp", nirrep_, nsopi_, nsopi_);

    // Set up DIIS
    DIISManager scfDiisManager(maxdiis_, "DCT DIIS Orbitals", DIISManager::LargestError, DIISManager::InCore);
    if ((nalpha_ + nbeta_) > 1) {
        scfDiisManager.set_error_vector_size(2, DIISEntry::Matrix, scf_error_a_.get(), DIISEntry::Matrix,
                                             scf_error_b_.get());
        scfDiisManager.set_vector_size(2, DIISEntry::Matrix, Fa_.get(), DIISEntry::Matrix, Fb_.get());
    }
    // Update the orbitals
    int nSCFCycles = 0;
    // Reset the booleans that control the convergence
    densityConverged_ = false;
    energyConverged_ = false;
    outfile->Printf("\tOrbital Updates\n");
    while ((!densityConverged_ || !orbitalsDone_ || !energyConverged_) && (nSCFCycles++ < maxiter_)) {
        std::string diisString;
        // Copy core hamiltonian into the Fock matrix array: F = H
        Fa_->copy(so_h_);
        Fb_->copy(so_h_);
        // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
        process_so_ints();
        // Save F0 = H + G * Kappa for the Fock intermediate update in lambda iterations
        moF0a_->copy(Fa_);
        moF0b_->copy(Fb_);
        moF0a_->transform(Ca_);
        moF0b_->transform(Cb_);
        // Save old SCF energy
        old_total_energy_ = new_total_energy_;
        // Add non-idempotent density contribution (Tau) to the Fock matrix: F += Gbar * Tau
        Fa_->add(g_tau_a_);
        Fb_->add(g_tau_b_);
        // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
        moFa_->copy(Fa_);
        moFb_->copy(Fb_);
        // Compute new SCF energy
        compute_scf_energy();
        // Check SCF convergence
        orbitals_convergence_ = compute_scf_error_vector();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        if (orbitals_convergence_ < diis_start_thresh_ && (nalpha_ + nbeta_) > 1) {
            if (scfDiisManager.add_entry(4, scf_error_a_.get(), scf_error_b_.get(), Fa_.get(), Fb_.get()))
                diisString += "S";
            if (scfDiisManager.subspace_size() > mindiisvecs_ && (nalpha_ + nbeta_) > 1) {
                diisString += "/E";
                scfDiisManager.extrapolate(2, Fa_.get(), Fb_.get());
            }
        }
        // Transform the Fock matrix to the symmetrically orhogonalized basis set and digonalize it
        // Obtain new orbitals
        Fa_->transform(s_half_inv_);
        Fa_->diagonalize(tmp, epsilon_a_);
        old_ca_->copy(Ca_);
        Ca_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        Fb_->transform(s_half_inv_);
        Fb_->diagonalize(tmp, epsilon_b_);
        old_cb_->copy(Cb_);
        Cb_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        // Make sure that the orbital phase is retained
        correct_mo_phases(false);
        // Update SCF density (Kappa) and check its RMS
        densityConverged_ = update_scf_density() < orbitals_threshold_;
        // Compute the DCT energy
        new_total_energy_ = scf_energy_ + lambda_energy_;
        // Check convergence of the total DCT energy
        energyConverged_ = std::fabs(new_total_energy_ - old_total_energy_) < energy_threshold_;
        outfile->Printf("\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n", nSCFCycles, orbitals_convergence_,
                        cumulant_convergence_, new_total_energy_ - old_total_energy_, new_total_energy_,
                        diisString.c_str());
        if (std::fabs(orbitals_convergence_) > 100.0) throw PSIEXCEPTION("DCT orbital updates diverged");
    }
    orbitalsDone_ = nSCFCycles == 1;
    energyConverged_ = false;
    // Transform the Fock matrix to the MO basis
    moFa_->transform(Ca_);
    moFb_->transform(Cb_);
    // Transform two-electron integrals to the MO basis using new orbitals, build denominators
    transform_integrals();
}

void DCTSolver::run_simult_dct() {
    // This is the simultaneous orbital/lambda update algorithm
    int cycle = 0;
    outfile->Printf(
        "\n\n\t*=================================================================================*\n"
        "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
        "\t*---------------------------------------------------------------------------------*\n");

    auto tmp = std::make_shared<Matrix>("temp", nirrep_, nsopi_, nsopi_);
    // Set up the DIIS manager
    DIISManager diisManager(maxdiis_, "DCT DIIS vectors");
    dpdbuf4 Laa, Lab, Lbb;
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Lambda <OO|VV>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Lambda <oo|vv>");
    diisManager.set_error_vector_size(5, DIISEntry::Matrix, scf_error_a_.get(), DIISEntry::Matrix, scf_error_b_.get(),
                                      DIISEntry::DPDBuf4, &Laa, DIISEntry::DPDBuf4, &Lab, DIISEntry::DPDBuf4, &Lbb);
    diisManager.set_vector_size(5, DIISEntry::Matrix, Fa_.get(), DIISEntry::Matrix, Fb_.get(), DIISEntry::DPDBuf4, &Laa,
                                DIISEntry::DPDBuf4, &Lab, DIISEntry::DPDBuf4, &Lbb);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);
    while ((!orbitalsDone_ || !cumulantDone_ || !densityConverged_ || !energyConverged_) && cycle++ < maxiter_) {
        std::string diisString;
        // Save the old energy
        old_total_energy_ = new_total_energy_;
        // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
        build_tau();
        if (exact_tau_) {
            refine_tau();
        }
        transform_tau();
        if (options_.get_str("DCT_TYPE") == "DF" && options_.get_str("AO_BASIS") == "NONE") {
            build_DF_tensors_UHF();

            auto mo_h_A = std::make_shared<Matrix>("MO-based H Alpha", nirrep_, nmopi_, nmopi_);
            mo_h_A->copy(so_h_);
            mo_h_A->transform(Ca_);

            auto mo_h_B = std::make_shared<Matrix>("MO-based H Beta", nirrep_, nmopi_, nmopi_);
            mo_h_B->copy(so_h_);
            mo_h_B->transform(Cb_);

            moFa_->copy(mo_h_A);
            moFb_->copy(mo_h_B);

            moFa_->add(mo_gbarGamma_A_);
            moFb_->add(mo_gbarGamma_B_);

            // Back-transform the Fock matrix to the SO basis: F_so = (Ct)^-1 F_mo C^-1 = (C^-1)t F_mo C^-1
            auto Ca_inverse = std::make_shared<Matrix>("Ca_ inverse", nirrep_, nmopi_, nsopi_);
            Ca_inverse->copy(Ca_);
            Ca_inverse->general_invert();
            Fa_->copy(moFa_);
            Fa_->transform(Ca_inverse);

            auto Cb_inverse = std::make_shared<Matrix>("Cb_ inverse", nirrep_, nmopi_, nsopi_);
            Cb_inverse->copy(Cb_);
            Cb_inverse->general_invert();
            Fb_->copy(moFa_);
            Fb_->transform(Cb_inverse);

        } else {
            // Copy core hamiltonian into the Fock matrix array: F = H
            Fa_->copy(so_h_);
            Fb_->copy(so_h_);
            // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
            process_so_ints();
            // Add non-idempotent density contribution (Tau) to the Fock matrix: F += Gbar * Tau
            Fa_->add(g_tau_a_);
            Fb_->add(g_tau_b_);
            // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
            moFa_->copy(Fa_);
            moFb_->copy(Fb_);
            // Transform the Fock matrix to the MO basis
            moFa_->transform(Ca_);
            moFb_->transform(Cb_);
        }

        // Compute new SCF energy
        compute_scf_energy();
        // Add SCF energy contribution to the total DCT energy
        new_total_energy_ = scf_energy_;
        // Check SCF convergence
        orbitals_convergence_ = compute_scf_error_vector();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        // Build G and F intermediates needed for the density cumulant residual equations and DCT energy computation
        build_cumulant_intermediates();
        // Compute the residuals for density cumulant equations
        cumulant_convergence_ = compute_cumulant_residual();
        if (std::fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCT density cumulant equations diverged");
        // Check convergence for density cumulant iterations
        cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
        // Update density cumulant tensor
        update_cumulant_jacobi();
        // Compute new DCT energy (lambda contribution)
        compute_dct_energy();
        // Add lambda energy to the DCT total energy
        new_total_energy_ += lambda_energy_;
        // Check convergence of the total DCT energy
        energyConverged_ = std::fabs(old_total_energy_ - new_total_energy_) < energy_threshold_;
        if (orbitals_convergence_ < diis_start_thresh_ && cumulant_convergence_ < diis_start_thresh_) {
            // Store the DIIS vectors
            dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
            global_dpd_->buf4_init(&Raa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                                   "R <OO|VV>");
            global_dpd_->buf4_init(&Rab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                                   "R <Oo|Vv>");
            global_dpd_->buf4_init(&Rbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                                   "R <oo|vv>");
            global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                                   "Lambda <OO|VV>");
            global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                                   "Lambda <Oo|Vv>");
            global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                                   "Lambda <oo|vv>");
            if (diisManager.add_entry(10, scf_error_a_.get(), scf_error_b_.get(), &Raa, &Rab, &Rbb, Fa_.get(),
                                      Fb_.get(), &Laa, &Lab, &Lbb)) {
                diisString += "S";
            }
            if (diisManager.subspace_size() > mindiisvecs_) {
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
        Fb_->transform(s_half_inv_);
        Fb_->diagonalize(tmp, epsilon_b_);
        old_cb_->copy(Cb_);
        Cb_->gemm(false, false, 1.0, s_half_inv_, tmp, 0.0);
        // Make sure that the orbital phase is retained
        if (!correct_mo_phases(false)) {
            outfile->Printf(
                "\t\tThere was a problem correcting the MO phases.\n"
                "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
        }
        // Transform two-electron integrals to the MO basis using new orbitals, build denominators
        transform_integrals();
        // Update SCF density (Kappa) and check its RMS
        densityConverged_ = update_scf_density() < orbitals_threshold_;
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
