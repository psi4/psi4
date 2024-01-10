/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"

#include <algorithm>
#include <cmath>

namespace psi {
namespace dct {

void DCTSolver::run_simult_dct_oo() {
    if (options_.get_bool("ODC_GUESS")) run_simult_dc_guess();

    // This is the simultaneous orbital/lambda update algorithm for the orbital-optimized methods
    int cycle = 0;

    outfile->Printf(
        "\n\n\t*=================================================================================*\n"
        "\t* Cycle   Max Orb Grad    RMS Lambda Error   delta E        Total Energy     DIIS *\n"
        "\t*---------------------------------------------------------------------------------*\n");

    // Copy the reference orbitals and to use them as the reference for the orbital rotation
    old_ca_->copy(Ca_);
    old_cb_->copy(Cb_);

    // Set up the DIIS manager
    DIISManager diisManager(maxdiis_, "DCT DIIS vectors");
    dpdbuf4 Laa, Lab, Lbb;
    global_dpd_->buf4_init(&Laa, PSIF_DCT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    diisManager.set_error_vector_size(orbital_gradient_a_.get(), orbital_gradient_b_.get(), &Laa, &Lab, &Lbb);
    diisManager.set_vector_size(Xtotal_a_.get(), Xtotal_b_.get(), &Laa, &Lab, &Lbb);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);

    while ((!orbitalsDone_ || !cumulantDone_ || !energyConverged_) && cycle++ < maxiter_) {
        std::string diisString;
        compute_SO_tau_U();

        if (options_.get_str("DCT_TYPE") == "DF" && options_.get_str("AO_BASIS") == "NONE") {
            build_DF_tensors_UHF();

            auto mo_h_A = Matrix("MO-based H Alpha", nirrep_, nmopi_, nmopi_);
            mo_h_A.copy(so_h_);
            mo_h_A.transform(Ca_);

            auto mo_h_B = Matrix("MO-based H Beta", nirrep_, nmopi_, nmopi_);
            mo_h_B.copy(so_h_);
            mo_h_B.transform(Cb_);

            moFa_->copy(mo_h_A);
            moFb_->copy(mo_h_B);

            moFa_->add(mo_gbarGamma_A_);
            moFb_->add(mo_gbarGamma_B_);
        } else {
            // Copy core hamiltonian into the Fock matrix array: F = H
            Fa_->copy(so_h_);
            Fb_->copy(so_h_);
            // Build the new Fock matrix from the SO integrals
            process_so_ints();
            // Back up the SO basis Fock before it is symmetrically orthogonalized to transform it to the MO basis
            moFa_->copy(Fa_);
            moFb_->copy(Fb_);
            // Transform the Fock matrix to the MO basis
            moFa_->transform(Ca_);
            moFb_->transform(Cb_);
        }
        // Compute new SCF energy
        compute_scf_energy();
        // Save the old energy
        old_total_energy_ = new_total_energy_;
        // Add SCF energy contribution to the total DCT energy
        new_total_energy_ = scf_energy_;
        // Build cumulant intermediates needed for the density cumulant residual equations and DCT energy computation
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
        // Compute orbital gradient and check convergence
        orbitals_convergence_ = compute_orbital_residual();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        // Check convergence of the total DCT energy
        energyConverged_ = std::fabs(old_total_energy_ - new_total_energy_) < energy_threshold_;
        // Compute the orbital rotation step using Jacobi method
        compute_orbital_rotation_jacobi();
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
                                   "Amplitude <OO|VV>");
            global_dpd_->buf4_init(&Lab, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                                   "Amplitude <Oo|Vv>");
            global_dpd_->buf4_init(&Lbb, PSIF_DCT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                                   "Amplitude <oo|vv>");
            if (diisManager.add_entry(orbital_gradient_a_.get(), orbital_gradient_b_.get(), &Raa, &Rab, &Rbb,
                                      Xtotal_a_.get(), Xtotal_b_.get(), &Laa, &Lab, &Lbb)) {
                diisString += "S";
            }

            if (diisManager.subspace_size() > mindiisvecs_) {
                diisString += "/E";
                diisManager.extrapolate(Xtotal_a_.get(), Xtotal_b_.get(), &Laa, &Lab, &Lbb);
            }
            global_dpd_->buf4_close(&Raa);
            global_dpd_->buf4_close(&Rab);
            global_dpd_->buf4_close(&Rbb);
            global_dpd_->buf4_close(&Laa);
            global_dpd_->buf4_close(&Lab);
            global_dpd_->buf4_close(&Lbb);
        }
        // Obtain new orbitals
        rotate_orbitals();
        // Make sure that the orbital phase is retained
        if (!correct_mo_phases(false)) {
            outfile->Printf(
                "\t\tThere was a problem correcting the MO phases.\n"
                "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
        }
        // Transform two-electron integrals to the MO basis using new orbitals, build denominators
        transform_integrals();
        // Update SCF density (Kappa)
        update_scf_density();
        // If we've performed enough lambda updates since the last orbitals
        // update, reset the counter so another SCF update is performed
        outfile->Printf("\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n", cycle, orbitals_convergence_,
                        cumulant_convergence_, new_total_energy_ - old_total_energy_, new_total_energy_,
                        diisString.c_str());
    }

    outfile->Printf("\t*=================================================================================*\n");
}

void DCTSolver::run_simult_dc_guess() {
    double lambda_conv = cumulant_threshold_;
    double orbital_conv = orbitals_threshold_;
    double energy_conv = energy_threshold_;

    cumulant_threshold_ = options_.get_double("GUESS_R_CONVERGENCE");
    orbitals_threshold_ = options_.get_double("GUESS_R_CONVERGENCE");
    energy_threshold_ = options_.get_double("GUESS_R_CONVERGENCE");
    orbital_optimized_ = false;

    std::string prefix = options_.get_str("DCT_TYPE") == "DF" ? "DF-" : "";
    outfile->Printf("\n\n\tComputing the guess using the %s%s functional", prefix.c_str(),
                    exact_tau_ ? "DC-12" : "DC-06");
    outfile->Printf("\n\tGuess energy, orbitals and cumulants will be converged to %4.3e",
                    options_.get_double("GUESS_R_CONVERGENCE"));
    if (options_.get_str("REFERENCE") == "RHF")
        run_simult_dct_RHF();
    else
        run_simult_dct();

    orbital_optimized_ = true;
    cumulantDone_ = false;
    orbitalsDone_ = false;

    cumulant_threshold_ = lambda_conv;
    orbitals_threshold_ = orbital_conv;
    energy_threshold_ = energy_conv;

    outfile->Printf("\n\tNow running the %s%s computation...", prefix.c_str(),
                    options_.get_str("DCT_FUNCTIONAL").c_str());
}

double DCTSolver::compute_orbital_residual() {
    dct_timer_on("DCTSolver::compute_orbital_residual()");

    dpdfile2 Xai, Xia;

    auto is_df = options_.get_str("DCT_TYPE") == "DF";

    // Compute the unrelaxed densities for the orbital gradient
    // For the conventional integral case, we want the entire density.
    // For the DF case, we want only the RIFIT terms, which happens to be the cumulant.
    compute_unrelaxed_density_OOOO(is_df);
    compute_unrelaxed_density_OOVV(is_df);
    compute_unrelaxed_density_OVOV(is_df);

    // Compute the OV part of the orbital gradient
    compute_orbital_gradient_OV(is_df);

    // Compute the VO part of the orbital gradient
    compute_orbital_gradient_VO(is_df);

    global_dpd_->file2_init(&Xia, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->file2_init(&Xai, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");

    auto temp = Matrix(&Xia);
    temp.subtract(Matrix(&Xai).transpose());
    temp.scale(2);
    double maxGradient = temp.absmax();
    orbital_gradient_a_->set_block(slices_.at("ACTIVE_OCC_A"), slices_.at("ACTIVE_VIR_A"), temp);
    temp.scale(-1);
    orbital_gradient_a_->set_block(slices_.at("ACTIVE_VIR_A"), slices_.at("ACTIVE_OCC_A"), temp.transpose());

    global_dpd_->file2_close(&Xai);
    global_dpd_->file2_close(&Xia);

    global_dpd_->file2_init(&Xia, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->file2_init(&Xai, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");

    temp = Matrix(&Xia);
    temp.subtract(Matrix(&Xai).transpose());
    temp.scale(2);
    maxGradient = std::max(maxGradient, temp.absmax());
    orbital_gradient_b_->set_block(slices_.at("ACTIVE_OCC_B"), slices_.at("ACTIVE_VIR_B"), temp);
    temp.scale(-1);
    orbital_gradient_b_->set_block(slices_.at("ACTIVE_VIR_B"), slices_.at("ACTIVE_OCC_B"), temp.transpose());

    global_dpd_->file2_close(&Xai);
    global_dpd_->file2_close(&Xia);

    dct_timer_off("DCTSolver::compute_orbital_residual()");

    return maxGradient;
}

void DCTSolver::compute_orbital_gradient_OV(bool separate_gbargamma) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdbuf4 L, W, LL;
    dpdfile2 X, H, T;
    dpdfile2 T_VV, T_vv;
    dpdfile2 Y2_OV, Y2_ov;

    // X_OV: One-electron contributions

    // X_IA = H_IB Tau_BA
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    global_dpd_->file2_init(&T, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->contract222(&H, &T, &X, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T);
    global_dpd_->file2_close(&H);
    global_dpd_->file2_close(&X);

    // X_ia = H_ib Tau_ba
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    global_dpd_->file2_init(&T, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
    global_dpd_->contract222(&H, &T, &X, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T);
    global_dpd_->file2_close(&H);
    global_dpd_->file2_close(&X);

    // X_OV: Two-electron contributions

    //
    // 2 * <OV||VV> Г_VVVV
    //
    // The VVVV block is expensive, so we avoid computing it by changing the contraction order.

    if (!separate_gbargamma) {
        // Compute contributions from VVVV density
        // 1. X_ia <-- <ib||cd> tau_ca tau_db
        global_dpd_->file2_init(&T_VV, PSIF_DCT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
        global_dpd_->file2_init(&T_vv, PSIF_DCT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        // Alpha contribution X_IA
        global_dpd_->file2_init(&Y2_OV, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "Y2 <O|V>");

        // Y2_IA = <IA|CD> Tau_CD
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 0,
                               "MO Ints <OV|VV>");
        global_dpd_->contract422(&I, &T_VV, &Y2_OV, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&I);
        // Y2_IA -= (IA|CD) Tau_CD
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V>=V]+"), 0,
                               "MO Ints (OV|VV)");
        global_dpd_->contract422(&I, &T_VV, &Y2_OV, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&I);
        // Y2_IA -= (IA|cd) Tau_cd
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"), ID("[O,V]"), ID("[v>=v]+"), 0,
                               "MO Ints (OV|vv)");
        global_dpd_->contract422(&I, &T_vv, &Y2_OV, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&I);

        // X_IA -= Y2_IC Tau_CA
        global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
        global_dpd_->contract222(&Y2_OV, &T_VV, &X, 0, 1, -1.0, 1.0);
        global_dpd_->file2_close(&X);

        global_dpd_->file2_close(&Y2_OV);

        // Beta contribution X_ia
        global_dpd_->file2_init(&Y2_ov, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "Y2 <o|v>");

        // Y2_ib = <ib|cd> Tau_cd
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"), ID("[o,v]"), ID("[v,v]"), 0,
                               "MO Ints <ov|vv>");
        global_dpd_->contract422(&I, &T_vv, &Y2_ov, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&I);
        // Y2_ib -= (ib|cd) Tau_cd
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"), ID("[o,v]"), ID("[v>=v]+"), 0,
                               "MO Ints (ov|vv)");
        global_dpd_->contract422(&I, &T_vv, &Y2_ov, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&I);
        // Y2_ib -= (ib|CD) Tau_CD
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[V>=V]+"), ID("[o,v]"), ID("[V>=V]+"), 0,
                               "MO Ints (ov|VV)");
        global_dpd_->contract422(&I, &T_VV, &Y2_ov, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&I);

        // X_ia -= Y2_ic Tau_ca
        global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
        global_dpd_->contract222(&Y2_ov, &T_vv, &X, 0, 1, -1.0, 1.0);
        global_dpd_->file2_close(&X);

        global_dpd_->file2_close(&Y2_ov);

        global_dpd_->file2_close(&T_VV);
        global_dpd_->file2_close(&T_vv);
    } else {
        // We already computed gbar gamma! We can use that intermediate now to compute ALL terms involving the
        // product of the 1RDM in our orbital gradient component. For this reason, we need to use the cumulant
        // and not the full gamma for the other 2RDM terms.
        // All we need is (gbargamma)^i_p gamma^p_a.
        auto zero = Dimension(nirrep_);
        const auto& SO = slices_.at("SO");
        const auto& MO = slices_.at("MO");
        // TODO: These slices grab blocks from the entire MO space.
        // When I have frozen core, mo_gammaA_ should definitely be active indices only.
        // This means I'll need a slice that grabs occ/vir from the active subspace, not the entire MO space.
        // I don't know what mo_gbarGamma_A_ should look like then.
        // I'll worry about how to implement this in the near future "add frozen core" PR.  JPM 04/23/21
        auto gbar_alpha_block = mo_gbarGamma_A_.get_block(slices_.at("ACTIVE_OCC_A"), MO);
        auto gbar_beta_block = mo_gbarGamma_B_.get_block(slices_.at("ACTIVE_OCC_B"), MO);
        auto gamma_alpha_block = mo_gammaA_.get_block(MO, slices_.at("ACTIVE_VIR_A"));
        auto gamma_beta_block = mo_gammaB_.get_block(MO, slices_.at("ACTIVE_VIR_B"));
        auto alpha_jk = linalg::doublet(gbar_alpha_block, gamma_alpha_block, false, false);
        auto beta_jk = linalg::doublet(gbar_beta_block, gamma_beta_block, false, false);

        global_dpd_->file2_init(&H, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X JK <O|V>");
        alpha_jk->write_to_dpdfile2(&H);
        global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
        global_dpd_->file2_axpy(&H, &X, 1.0, 0);
        global_dpd_->file2_close(&X);
        global_dpd_->file2_close(&H);

        global_dpd_->file2_init(&H, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X JK <o|v>");
        beta_jk->write_to_dpdfile2(&H);
        global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
        global_dpd_->file2_axpy(&H, &X, 1.0, 0);
        global_dpd_->file2_close(&X);
        global_dpd_->file2_close(&H);
    }

    // 2. X_ia <-- 1/4 <ib||cd> lambda_abkl lambda_klcd

    //  W_IBKL = <IB||CD> lambda_KLCD
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V>V]-"), ID("[O,V]"), ID("[V,V]"), 1,
                           "MO Ints <OV|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V>V]-"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,O]"), ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    global_dpd_->contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&W);

    //  W_KlIb = 2 lambda_KlCd <Ib|Cd>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"), ID("[O,v]"), ID("[V,v]"), 0,
                           "MO Ints <Ov|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,v]"), ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    global_dpd_->contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&W);

    //  W_LkBi = 2 lambda_LkDc <Bi|Dc>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"), ID("[V,o]"), ID("[V,v]"), 0,
                           "MO Ints <Vo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,o]"), ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    global_dpd_->contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&W);

    //  W_ibkl = <ib||cd> lambda_klcd
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v>v]-"), ID("[o,v]"), ID("[v,v]"), 1,
                           "MO Ints <ov|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v>v]-"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,o]"), ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    global_dpd_->contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&W);

    // X_IA +=  1/4 W_IBKL lambda_KLAB
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,V]"), ID("[O,O]"), ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    global_dpd_->buf4_init(&LL, PSIF_DCT_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           "Amplitude <OO|VV>");

    global_dpd_->contract442(&W, &LL, &X, 0, 2, 0.25, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LL);
    global_dpd_->file2_close(&X);

    // X_IA +=  1/2 W_KlIb lambda_KlAb
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[O,v]"), ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    global_dpd_->buf4_init(&LL, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");

    global_dpd_->contract442(&W, &LL, &X, 2, 2, 0.5, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LL);
    global_dpd_->file2_close(&X);

    // X_ia +=  1/2 W_LkBi lambda_LkBa
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,o]"), ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    global_dpd_->buf4_init(&LL, PSIF_DCT_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           "Amplitude <Oo|Vv>");

    global_dpd_->contract442(&W, &LL, &X, 3, 3, 0.5, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LL);
    global_dpd_->file2_close(&X);

    // X_ia += 1/4 W_ibkl lambda_klab
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&W, PSIF_DCT_DPD, 0, ID("[o,v]"), ID("[o,o]"), ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    global_dpd_->buf4_init(&LL, PSIF_DCT_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           "Amplitude <oo|vv>");

    global_dpd_->contract442(&W, &LL, &X, 0, 2, 0.25, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LL);
    global_dpd_->file2_close(&X);

    //
    // <OO||OV> Г_OOVV
    //

    std::string density_variable = separate_gbargamma ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    // X_IA += <BI||JK> Г_BAJK
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), 1,
                           "MO Ints <VO|OO>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"), ID("[V>V]-"), ID("[O>O]-"), 0,
                           varname("<VV|OO>"));

    global_dpd_->contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_IA += 2 * <Ib|Jk> Г_AbJk
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"), ID("[O,v]"), ID("[O,o]"), 0,
                           "MO Ints <Ov|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), 0,
                           varname("<Vv|Oo>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia += <bi||jk> Г_bajk
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"), ID("[v,o]"), ID("[o,o]"), 1,
                           "MO Ints <vo|oo>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[v,v]"), ID("[o,o]"), ID("[v>v]-"), ID("[o>o]-"), 0,
                           varname("<vv|oo>"));

    global_dpd_->contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia += 2 * <Bi|Jk> Г_BaJk
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"), ID("[V,o]"), ID("[O,o]"), 0,
                           "MO Ints <Vo|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), 0,
                           varname("<Vv|Oo>"));

    global_dpd_->contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <OO||OV> Г_OVOV
    //

    // X_IA += <JB||KI> Г_JBKA
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"), ID("[O,V]"), ID("[O,O]"), 1,
                           "MO Ints <OV|OO>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("<OV|OV>"));

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_IA += <kB|jI> Г_kBjA
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"), ID("[o,V]"), ID("[o,O]"), 0,
                           "MO Ints <oV|oO>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0,
                           varname("<oV|oV>"));

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_IA -= <Kb|jI> Г_KbjA
    // Note: <Kb|jI> integrals are resorted <bK|jI> integrals.
    // <Kb||jI> Г_KbjA = (-1) * <bK|jI> * Г_KbjA
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,O]"), ID("[O,v]"), ID("[o,O]"), 0,
                           "MO Ints <Ov|oO>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0,
                           varname("<Ov|oV>"));

    global_dpd_->contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia += <jb||ki> Г_jbka
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"), ID("[o,v]"), ID("[o,o]"), 1,
                           "MO Ints <ov|oo>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           varname("<ov|ov>"));

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia += <Kb|Ji> Г_KbJa
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"), ID("[O,v]"), ID("[O,o]"), 0,
                           "MO Ints <Ov|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0,
                           varname("<Ov|Ov>"));

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia -= <kB|Ji> Г_kBJa
    // Note: <kB|Ji> integrals are resorted <Bk|Ji> integrals.
    // <kB||Ji> Г_kBJa = (-1) * <Bk|Ji> * Г_kBJa

    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[O,o]"), ID("[o,V]"), ID("[O,o]"), 0,
                           "MO Ints <oV|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), 0,
                           varname("<oV|Ov>"));

    global_dpd_->contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    psio_->close(PSIF_DCT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void DCTSolver::compute_orbital_gradient_VO(bool separate_gbargamma) {
    psio_->open(PSIF_DCT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdfile2 X, H, T;

    // X_VO: One-electron contributions

    // X_AI = H_JA Tau_JI
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    global_dpd_->file2_init(&T, PSIF_DCT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_mat_init(&X);
    global_dpd_->file2_mat_init(&H);
    global_dpd_->file2_mat_init(&T);
    global_dpd_->file2_mat_rd(&H);
    global_dpd_->file2_mat_rd(&T);
    for (int h = 0; h < nirrep_; ++h) {
#pragma omp parallel for
        for (int i = 0; i < naoccpi_[h]; ++i) {
            for (int a = 0; a < navirpi_[h]; ++a) {
                double value = 0.0;
                for (int j = 0; j < naoccpi_[h]; ++j) {
                    value += H.matrix[h][j][a] * (T.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                }
                X.matrix[h][a][i] = value;
            }
        }
    }
    global_dpd_->file2_mat_wrt(&X);
    global_dpd_->file2_close(&T);
    global_dpd_->file2_close(&H);
    global_dpd_->file2_close(&X);

    // X_ai = H_ja Tau_ji
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    global_dpd_->file2_init(&T, PSIF_DCT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_mat_init(&X);
    global_dpd_->file2_mat_init(&H);
    global_dpd_->file2_mat_init(&T);
    global_dpd_->file2_mat_rd(&H);
    global_dpd_->file2_mat_rd(&T);
    for (int h = 0; h < nirrep_; ++h) {
#pragma omp parallel for
        for (int i = 0; i < nboccpi_[h]; ++i) {
            for (int a = 0; a < nbvirpi_[h]; ++a) {
                double value = 0.0;
                for (int j = 0; j < nboccpi_[h]; ++j) {
                    value += H.matrix[h][j][a] * (T.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                }
                X.matrix[h][a][i] = value;
            }
        }
    }
    global_dpd_->file2_mat_wrt(&X);
    global_dpd_->file2_close(&T);
    global_dpd_->file2_close(&H);
    global_dpd_->file2_close(&X);

    if (separate_gbargamma) {
        // We already computed gbar gamma! We can use that intermediate now to compute ALL terms involving the
        // product of the 1RDM in our orbital gradient component. For this reason, we need to use the cumulant
        // and not the full gamma for the other 2RDM terms.
        // All we need is (gbargamma)^a_p gamma^p_i.
        auto zero = Dimension(nirrep_);
        auto gbar_alpha_block = mo_gbarGamma_A_.get_block(Slice(nalphapi_, nmopi_), Slice(zero, nmopi_));
        auto gbar_beta_block = mo_gbarGamma_B_.get_block(Slice(nbetapi_, nmopi_), Slice(zero, nmopi_));
        auto gamma_alpha_block = mo_gammaA_.get_block(Slice(zero, nmopi_), Slice(zero, nalphapi_));
        auto gamma_beta_block = mo_gammaB_.get_block(Slice(zero, nmopi_), Slice(zero, nbetapi_));
        auto alpha_jk = linalg::doublet(gbar_alpha_block, gamma_alpha_block, false, false);
        auto beta_jk = linalg::doublet(gbar_beta_block, gamma_beta_block, false, false);

        global_dpd_->file2_init(&H, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X JK <V|O>");
        alpha_jk->write_to_dpdfile2(&H);
        global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
        global_dpd_->file2_axpy(&H, &X, 1.0, 0);
        global_dpd_->file2_close(&X);
        global_dpd_->file2_close(&H);

        global_dpd_->file2_init(&H, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X JK <v|o>");
        beta_jk->write_to_dpdfile2(&H);
        global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
        global_dpd_->file2_axpy(&H, &X, 1.0, 0);
        global_dpd_->file2_close(&X);
        global_dpd_->file2_close(&H);
    }

    // X_VO: Two-electron contributions

    //
    // 2 * <VO||OO> Г_OOOO
    //

    std::string density_variable = separate_gbargamma ? "Lambda " : "Gamma ";
    auto varname = [&density_variable](const std::string& x) { return (density_variable + x); };

    // X_AI += 2 * <AJ||KL> Г_IJKL
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"), ID("[V,O]"), ID("[O,O]"), 1,
                           "MO Ints <VO|OO>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"), ID("[O>O]-"), ID("[O>O]-"), 0,
                           varname("<OO|OO>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI += 4 * <Aj|Kl> Г_IjKl
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"), ID("[V,o]"), ID("[O,o]"), 0,
                           "MO Ints <Vo|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0,
                           varname("<Oo|Oo>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += 2 * <aj||kl> Г_ijkl
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"), ID("[v,o]"), ID("[o,o]"), 1,
                           "MO Ints <vo|oo>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"), ID("[o>o]-"), ID("[o>o]-"), 0,
                           varname("<oo|oo>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI += 4 * <Ja|Kl> Г_JiKl
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"), ID("[O,v]"), ID("[O,o]"), 0,
                           "MO Ints <Ov|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), ID("[O,o]"), 0,
                           varname("<Oo|Oo>"));

    global_dpd_->contract442(&I, &G, &X, 1, 1, 4.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <VO||VV> Г_OOVV
    //

    // X_AI += <JA||BC> Г_JIBC
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 1,
                           "MO Ints <OV|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>O]-"), ID("[V>V]-"), 0,
                           varname("<OO|VV>"));

    global_dpd_->contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI += <Aj|Bc> Г_IjBc
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"), ID("[V,o]"), ID("[V,v]"), 0,
                           "MO Ints <Vo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           varname("<Oo|Vv>"));

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += <ja||bc> Г_jibc
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"), ID("[o,v]"), ID("[v,v]"), 1,
                           "MO Ints <ov|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"), ID("[o>o]-"), ID("[v>v]-"), 0,
                           varname("<oo|vv>"));

    global_dpd_->contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += <Ja|Bc> Г_JiBc
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"), ID("[O,v]"), ID("[V,v]"), 0,
                           "MO Ints <Ov|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                           varname("<Oo|Vv>"));

    global_dpd_->contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <OV||VV> Г_OVOV
    //

    // X_AI += <JB||AC> Г_JBIC
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V,V]"), 1,
                           "MO Ints <OV|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0,
                           varname("<OV|OV>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI += <Jb|Ac> Г_JbIc
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"), ID("[O,v]"), ID("[V,v]"), 0,
                           "MO Ints <Ov|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), ID("[O,v]"), 0,
                           varname("<Ov|Ov>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI -= <jB|Ac> Г_jBIc
    // Note: <jB|Ac> integrals are resorted <Bj|Ac> integrals.
    // <jB||Ac> Г_jBIc = (-1) * <Bj|Ac> Г_jBIc
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[V,v]"), ID("[o,V]"), ID("[V,v]"), 0,
                           "MO Ints <oV|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), 0,
                           varname("<oV|Ov>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += <jb||ac> Г_jbic
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"), ID("[o,v]"), ID("[v,v]"), 1,
                           "MO Ints <ov|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), ID("[o,v]"), 0,
                           varname("<ov|ov>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += <jB|aC> Г_jBiC
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"), ID("[o,V]"), ID("[v,V]"), 0,
                           "MO Ints <oV|vV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), ID("[o,V]"), 0,
                           varname("<oV|oV>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai -= <Jb|aC> Г_JbiC
    // Note: <Jb|aC> integrals are resorted <bJ|aC> integrals.
    // <Jb||aC> Г_JbiC = (-1) * <bJ|aC> Г_JbiC
    global_dpd_->file2_init(&X, PSIF_DCT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[v,V]"), ID("[O,v]"), ID("[v,V]"), 0,
                           "MO Ints <Ov|vV>");
    global_dpd_->buf4_init(&G, PSIF_DCT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"), ID("[O,v]"), ID("[o,V]"), 0,
                           varname("<Ov|oV>"));

    global_dpd_->contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    psio_->close(PSIF_DCT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void DCTSolver::compute_orbital_rotation_jacobi() {
    dct_timer_on("DCTSolver::compute_orbital_rotation_jacobi()");

    // Determine the orbital rotation step
    // Alpha spin
    auto X_a = Matrix("Alpha orbital step", nirrep_, nmopi_, nmopi_);
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < naoccpi_[h]; ++i) {
            for (int a = naoccpi_[h]; a < nmopi_[h]; ++a) {
                double value = orbital_gradient_a_->get(h, i, a) /
                               (2.0 * (moFa_->get(h, i, i) - moFa_->get(h, a, a)) + orbital_level_shift_);
                X_a.set(h, i, a, value);
                X_a.set(h, a, i, (-1.0) * value);
            }
        }
    }

    // Beta spin
    auto X_b = Matrix("Beta orbital step", nirrep_, nmopi_, nmopi_);
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < nboccpi_[h]; ++i) {
            for (int a = nboccpi_[h]; a < nmopi_[h]; ++a) {
                double value = orbital_gradient_b_->get(h, i, a) /
                               (2.0 * (moFb_->get(h, i, i) - moFb_->get(h, a, a)) + orbital_level_shift_);
                X_b.set(h, i, a, value);
                X_b.set(h, a, i, (-1.0) * value);
            }
        }
    }

    // Determine the rotation generator with respect to the reference orbitals
    Xtotal_a_->add(X_a);
    Xtotal_b_->add(X_b);

    dct_timer_off("DCTSolver::compute_orbital_rotation_jacobi()");
}

void DCTSolver::rotate_matrix(const Matrix& X, const Matrix& old, Matrix& target) {
    auto U = *X.clone();
    U.expm(4, true);
    target.gemm(false, false, 1.0, old, U, 0.0);
}

void DCTSolver::rotate_orbitals() {
    dct_timer_on("DCTSolver::rotate_orbitals()");

    rotate_matrix(*Xtotal_a_, *old_ca_, *Ca_);
    rotate_matrix(*Xtotal_b_, *old_cb_, *Cb_);

    dct_timer_off("DCTSolver::rotate_orbitals()");
}

}  // namespace dct
}  // namespace psi
