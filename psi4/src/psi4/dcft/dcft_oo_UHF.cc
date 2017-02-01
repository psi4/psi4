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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdiis/diismanager.h"
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

void
DCFTSolver::run_simult_dcft_oo()
{

    if (options_.get_bool("ODC_GUESS")) run_simult_dc_guess();

    // This is the simultaneous orbital/lambda update algorithm for the orbital-optimized methods
    int cycle = 0;

    outfile->Printf( "\n\n\t*=================================================================================*\n"
                         "\t* Cycle   Max Orb Grad    RMS Lambda Error   delta E        Total Energy     DIIS *\n"
                         "\t*---------------------------------------------------------------------------------*\n");

    // Copy the reference orbitals and to use them as the reference for the orbital rotation
    old_ca_->copy(Ca_);
    old_cb_->copy(Cb_);

    // Set up the DIIS manager
    DIISManager diisManager(maxdiis_, "DCFT DIIS vectors");
    dpdbuf4 Laa, Lab, Lbb;
    global_dpd_->buf4_init(&Laa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&Lbb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    diisManager.set_error_vector_size(5, DIISEntry::Matrix, orbital_gradient_a_.get(),
                                         DIISEntry::Matrix, orbital_gradient_b_.get(),
                                         DIISEntry::DPDBuf4, &Laa,
                                         DIISEntry::DPDBuf4, &Lab,
                                         DIISEntry::DPDBuf4, &Lbb);
    diisManager.set_vector_size(5, DIISEntry::Matrix, Xtotal_a_.get(),
                                   DIISEntry::Matrix, Xtotal_b_.get(),
                                   DIISEntry::DPDBuf4, &Laa,
                                   DIISEntry::DPDBuf4, &Lab,
                                   DIISEntry::DPDBuf4, &Lbb);
    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Lbb);

    while((!orbitalsDone_ || !cumulantDone_ || !densityConverged_ || !energyConverged_)
            && cycle++ < maxiter_){
        std::string diisString;
        // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
        build_tau();
        if (exact_tau_) {
            refine_tau();
        }
        transform_tau();

        if (options_.get_str("DCFT_TYPE") == "DF" && options_.get_str("AO_BASIS") == "NONE"){

            build_DF_tensors_UHF();

            SharedMatrix mo_h_A = SharedMatrix(new Matrix("MO-based H Alpha", nirrep_, nmopi_, nmopi_));
            mo_h_A->copy(so_h_);
            mo_h_A->transform(Ca_);

            SharedMatrix mo_h_B = SharedMatrix(new Matrix("MO-based H Beta", nirrep_, nmopi_, nmopi_));
            mo_h_B->copy(so_h_);
            mo_h_B->transform(Cb_);

            moFa_->copy(mo_h_A);
            moFb_->copy(mo_h_B);

            moFa_->add(mo_gbarGamma_A_);
            moFb_->add(mo_gbarGamma_B_);
        }
        else{
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
        // Save the old energy
        old_total_energy_ = new_total_energy_;
        // Add SCF energy contribution to the total DCFT energy
        new_total_energy_ = scf_energy_;
        // Build cumulant intermediates needed for the density cumulant residual equations and DCFT energy computation
        build_cumulant_intermediates();
        // Compute the residuals for density cumulant equations
        cumulant_convergence_ = compute_cumulant_residual();
        if (fabs(cumulant_convergence_) > 100.0) throw PSIEXCEPTION("DCFT density cumulant equations diverged");
        // Check convergence for density cumulant iterations
        cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
        // Update density cumulant tensor
        update_cumulant_jacobi();
        // Compute new DCFT energy (lambda contribution)
        compute_dcft_energy();
        // Add lambda energy to the DCFT total energy
        new_total_energy_ += lambda_energy_;
        // Compute orbital gradient and check convergence
        orbitals_convergence_ = compute_orbital_residual();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        // Check convergence of the total DCFT energy
        energyConverged_ = fabs(old_total_energy_ - new_total_energy_) < energy_threshold_;
        // Compute the orbital rotation step using Jacobi method
        compute_orbital_rotation_jacobi();
        if(orbitals_convergence_ < diis_start_thresh_ && cumulant_convergence_ < diis_start_thresh_){
            //Store the DIIS vectors
            dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
            global_dpd_->buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
            global_dpd_->buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
            global_dpd_->buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
            global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
            global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
            global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
            if(diisManager.add_entry(10, orbital_gradient_a_.get(), orbital_gradient_b_.get(), &Raa, &Rab, &Rbb,
                                     Xtotal_a_.get(), Xtotal_b_.get(), &Laa, &Lab, &Lbb)){
                diisString += "S";
            }
            if(diisManager.subspace_size() > mindiisvecs_){
                diisString += "/E";
                diisManager.extrapolate(5, Xtotal_a_.get(), Xtotal_b_.get(), &Laa, &Lab, &Lbb);
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
        if(!correct_mo_phases(false)){
            outfile->Printf("\t\tThere was a problem correcting the MO phases.\n"
                            "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
        }
        // Transform two-electron integrals to the MO basis using new orbitals, build denominators
        transform_integrals();
        // Update SCF density (Kappa) and check its RMS
        densityConverged_ = update_scf_density() < orbitals_threshold_;
        // If we've performed enough lambda updates since the last orbitals
        // update, reset the counter so another SCF update is performed
        outfile->Printf( "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                cycle, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                new_total_energy_, diisString.c_str());

    }

    outfile->Printf( "\t*=================================================================================*\n");
}

void
DCFTSolver::run_simult_dc_guess()
{

    double lambda_conv = cumulant_threshold_;
    double orbital_conv = orbitals_threshold_;
    double energy_conv = energy_threshold_;

    cumulant_threshold_ = options_.get_double("GUESS_R_CONVERGENCE");
    orbitals_threshold_ = options_.get_double("GUESS_R_CONVERGENCE");
    energy_threshold_   = options_.get_double("GUESS_R_CONVERGENCE");
    orbital_optimized_ = false;

    std::string prefix = options_.get_str("DCFT_TYPE") == "DF"? "DF-" : "";
    outfile->Printf( "\n\n\tComputing the guess using the %s%s functional", prefix.c_str(), exact_tau_ ? "DC-12" : "DC-06");
    outfile->Printf( "\n\tGuess energy, orbitals and cumulants will be converged to %4.3e", options_.get_double("GUESS_R_CONVERGENCE"));
    if(options_.get_str("REFERENCE") == "RHF")
        run_simult_dcft_RHF();
    else
        run_simult_dcft();

    orbital_optimized_ = true;
    cumulantDone_ = false;
    orbitalsDone_ = false;

    cumulant_threshold_ = lambda_conv;
    orbitals_threshold_ = orbital_conv;
    energy_threshold_ = energy_conv;

    outfile->Printf( "\n\tNow running the %s%s computation...", prefix.c_str(), options_.get_str("DCFT_FUNCTIONAL").c_str());

}

double
DCFTSolver::compute_orbital_residual() {

    dcft_timer_on("DCFTSolver::compute_orbital_residual()");

    dpdfile2 Xai, Xia;

    // Compute the unrelaxed densities for the orbital gradient
    compute_unrelaxed_density_OOOO();
    compute_unrelaxed_density_OOVV();
    compute_unrelaxed_density_OVOV();

    // Compute the OV part of the orbital gradient
    compute_orbital_gradient_OV();

    // Compute the VO part of the orbital gradient
    compute_orbital_gradient_VO();

    global_dpd_->file2_init(&Xia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->file2_init(&Xai, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->file2_mat_init(&Xia);
    global_dpd_->file2_mat_init(&Xai);
    global_dpd_->file2_mat_rd(&Xia);
    global_dpd_->file2_mat_rd(&Xai);

    double maxGradient = 0.0;
    // Alpha spin
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                double value = 2.0 * (Xia.matrix[h][i][a] - Xai.matrix[h][a][i]);
                maxGradient = (fabs(value) > maxGradient) ? fabs(value) : maxGradient;
                orbital_gradient_a_->set(h, i, a + naoccpi_[h], value);
                orbital_gradient_a_->set(h, a + naoccpi_[h], i, (-1.0) * value);
            }
        }
    }

    global_dpd_->file2_close(&Xai);
    global_dpd_->file2_close(&Xia);

    global_dpd_->file2_init(&Xia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->file2_init(&Xai, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->file2_mat_init(&Xia);
    global_dpd_->file2_mat_init(&Xai);
    global_dpd_->file2_mat_rd(&Xia);
    global_dpd_->file2_mat_rd(&Xai);

    // Beta spin
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = 0; a < nbvirpi_[h]; ++a){
                double value = 2.0 * (Xia.matrix[h][i][a] - Xai.matrix[h][a][i]);
                maxGradient = (fabs(value) > maxGradient) ? fabs(value) : maxGradient;
                orbital_gradient_b_->set(h, i, a + nboccpi_[h], value);
                orbital_gradient_b_->set(h, a + nboccpi_[h], i, (-1.0) * value);
            }
        }
    }

    global_dpd_->file2_close(&Xai);
    global_dpd_->file2_close(&Xia);

    dcft_timer_off("DCFTSolver::compute_orbital_residual()");

    return maxGradient;
}

void
DCFTSolver::compute_orbital_gradient_OV() {

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdbuf4 L, W, LL;
    dpdfile2 X, H, T;
    dpdfile2 T_VV, T_vv;
    dpdfile2 Y2_OV, Y2_ov;

    // X_OV: One-electron contributions

    // X_IA = H_IB Tau_BA
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    global_dpd_->file2_init(&T, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->contract222(&H, &T, &X, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T);
    global_dpd_->file2_close(&H);
    global_dpd_->file2_close(&X);

    // X_ia = H_ib Tau_ba
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    global_dpd_->file2_init(&T, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
    global_dpd_->contract222(&H, &T, &X, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&T);
    global_dpd_->file2_close(&H);
    global_dpd_->file2_close(&X);

    // X_OV: Two-electron contributions

    //
    // 2 * <OV||VV> Г_VVVV
    //

    // Compute contributions from VVVV density
    // 1. X_ia <-- <ib||cd> tau_ca tau_db
    global_dpd_->file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    global_dpd_->file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    // Alpha contribution X_IA
    global_dpd_->file2_init(&Y2_OV, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Y2 <O|V>");

    // Y2_IA = <IA|CD> Tau_CD
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    global_dpd_->contract422(&I, &T_VV, &Y2_OV, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&I);
    // Y2_IA -= (IA|CD) Tau_CD
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    global_dpd_->contract422(&I, &T_VV, &Y2_OV, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // Y2_IA -= (IA|cd) Tau_cd
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    global_dpd_->contract422(&I, &T_vv, &Y2_OV, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // X_IA -= Y2_IC Tau_CA
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->contract222(&Y2_OV, &T_VV, &X, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&X);

    global_dpd_->file2_close(&Y2_OV);

    // Beta contribution X_ia
    global_dpd_->file2_init(&Y2_ov, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Y2 <o|v>");

    // Y2_ib = <ib|cd> Tau_cd
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
    global_dpd_->contract422(&I, &T_vv, &Y2_ov, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&I);
    // Y2_ib -= (ib|cd) Tau_cd
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
    global_dpd_->contract422(&I, &T_vv, &Y2_ov, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);
    // Y2_ib -= (ib|CD) Tau_CD
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[V>=V]+"),
                  ID("[o,v]"), ID("[V>=V]+"), 0, "MO Ints (ov|VV)");
    global_dpd_->contract422(&I, &T_VV, &Y2_ov, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&I);

    // X_ia -= Y2_ic Tau_ca
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->contract222(&Y2_ov, &T_vv, &X, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&X);

    global_dpd_->file2_close(&Y2_ov);

    global_dpd_->file2_close(&T_VV);
    global_dpd_->file2_close(&T_vv);

    // 2. X_ia <-- 1/4 <ib||cd> lambda_abkl lambda_klcd

    //  W_IBKL = <IB||CD> lambda_KLCD
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V>V]-"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    global_dpd_->contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&W);


    //  W_KlIb = 2 lambda_KlCd <Ib|Cd>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    global_dpd_->contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&W);

    //  W_LkBi = 2 lambda_LkDc <Bi|Dc>
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    global_dpd_->contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&W);

    //  W_ibkl = <ib||cd> lambda_klcd
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v>v]-"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    global_dpd_->contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&L);
    global_dpd_->buf4_close(&W);


    // X_IA +=  1/4 W_IBKL lambda_KLAB
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    global_dpd_->buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");

    global_dpd_->contract442(&W, &LL, &X, 0, 2, 0.25, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LL);
    global_dpd_->file2_close(&X);

    // X_IA +=  1/2 W_KlIb lambda_KlAb
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    global_dpd_->buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

    global_dpd_->contract442(&W, &LL, &X, 2, 2, 0.5, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LL);
    global_dpd_->file2_close(&X);

    // X_ia +=  1/2 W_LkBi lambda_LkBa
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    global_dpd_->buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

    global_dpd_->contract442(&W, &LL, &X, 3, 3, 0.5, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LL);
    global_dpd_->file2_close(&X);

    // X_ia += 1/4 W_ibkl lambda_klab
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    global_dpd_->buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");

    global_dpd_->contract442(&W, &LL, &X, 0, 2, 0.25, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&LL);
    global_dpd_->file2_close(&X);

    //
    // <OO||OV> Г_OOVV
    //

    // X_IA += <BI||JK> Г_BAJK
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V>V]-"), ID("[O>O]-"), 0, "Gamma <VV|OO>");

    global_dpd_->contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_IA += 2 * <Ib|Jk> Г_AbJk
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia += <bi||jk> Г_bajk
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v>v]-"), ID("[o>o]-"), 0, "Gamma <vv|oo>");

    global_dpd_->contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia += 2 * <Bi|Jk> Г_BaJk
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");

    global_dpd_->contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <OO||OV> Г_OVOV
    //

    // X_IA += <JB||KI> Г_JBKA
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 1, "MO Ints <OV|OO>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_IA += <kB|jI> Г_kBjA
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|oO>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_IA -= <Kb|jI> Г_KbjA
    // Note: <Kb|jI> integrals are resorted <bK|jI> integrals.
    // <Kb||jI> Г_KbjA = (-1) * <bK|jI> * Г_KbjA
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,O]"),
                  ID("[O,v]"), ID("[o,O]"), 0, "MO Ints <Ov|oO>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    global_dpd_->contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia += <jb||ki> Г_jbka
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 1, "MO Ints <ov|oo>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia += <Kb|Ji> Г_KbJa
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    global_dpd_->contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ia -= <kB|Ji> Г_kBJa
    // Note: <kB|Ji> integrals are resorted <Bk|Ji> integrals.
    // <kB||Ji> Г_kBJa = (-1) * <Bk|Ji> * Г_kBJa

    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[O,o]"),
                  ID("[o,V]"), ID("[O,o]"), 0, "MO Ints <oV|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "Gamma <oV|Ov>");

    global_dpd_->contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::compute_orbital_gradient_VO() {

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdfile2 X, H, T;

    // X_VO: One-electron contributions

    // X_AI = H_JA Tau_JI
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    global_dpd_->file2_init(&T, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    global_dpd_->file2_mat_init(&X);
    global_dpd_->file2_mat_init(&H);
    global_dpd_->file2_mat_init(&T);
    global_dpd_->file2_mat_rd(&H);
    global_dpd_->file2_mat_rd(&T);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < naoccpi_[h]; ++j){
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
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    global_dpd_->file2_init(&T, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    global_dpd_->file2_mat_init(&X);
    global_dpd_->file2_mat_init(&H);
    global_dpd_->file2_mat_init(&T);
    global_dpd_->file2_mat_rd(&H);
    global_dpd_->file2_mat_rd(&T);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < nboccpi_[h]; ++j){
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

    // X_VO: Two-electron contributions

    //
    // 2 * <VO||OO> Г_OOOO
    //

    // X_AI += 2 * <AJ||KL> Г_IJKL
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>O]-"), ID("[O>O]-"), 0, "Gamma <OO|OO>");

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI += 4 * <Aj|Kl> Г_IjKl
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    global_dpd_->contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += 2 * <aj||kl> Г_ijkl
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>o]-"), ID("[o>o]-"), 0, "Gamma <oo|oo>");

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI += 4 * <Ja|Kl> Г_JiKl
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    global_dpd_->contract442(&I, &G, &X, 1, 1, 4.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <VO||VV> Г_OOVV
    //

    // X_AI += <JA||BC> Г_JIBC
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Gamma <OO|VV>");

    global_dpd_->contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI += <Aj|Bc> Г_IjBc
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    global_dpd_->contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += <ja||bc> Г_jibc
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Gamma <oo|vv>");

    global_dpd_->contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += <Ja|Bc> Г_JiBc
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    global_dpd_->contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    //
    // <OV||VV> Г_OVOV
    //

    // X_AI += <JB||AC> Г_JBIC
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI += <Jb|Ac> Г_JbIc
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_AI -= <jB|Ac> Г_jBIc
    // Note: <jB|Ac> integrals are resorted <Bj|Ac> integrals.
    // <jB||Ac> Г_jBIc = (-1) * <Bj|Ac> Г_jBIc
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[V,v]"),
                  ID("[o,V]"), ID("[V,v]"), 0, "MO Ints <oV|Vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "Gamma <oV|Ov>");

    global_dpd_->contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += <jb||ac> Г_jbic
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai += <jB|aC> Г_jBiC
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    global_dpd_->contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    // X_ai -= <Jb|aC> Г_JbiC
    // Note: <Jb|aC> integrals are resorted <bJ|aC> integrals.
    // <Jb||aC> Г_JbiC = (-1) * <bJ|aC> Г_JbiC
    global_dpd_->file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[v,V]"),
                  ID("[O,v]"), ID("[v,V]"), 0, "MO Ints <Ov|vV>");
    global_dpd_->buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    global_dpd_->contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_close(&I);
    global_dpd_->file2_close(&X);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void
DCFTSolver::compute_orbital_rotation_jacobi() {

    dcft_timer_on("DCFTSolver::compute_orbital_rotation_jacobi()");

    // Determine the orbital rotation step
    // Alpha spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = naoccpi_[h]; a < nmopi_[h]; ++a){
                double value = orbital_gradient_a_->get(h, i, a) / (2.0 * (moFa_->get(h, i, i) - moFa_->get(h, a, a)) + orbital_level_shift_);
                X_a_->set(h, i, a, value);
                X_a_->set(h, a, i, (-1.0) * value);
            }
        }
    }

    // Beta spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = nboccpi_[h]; a < nmopi_[h]; ++a){
                double value = orbital_gradient_b_->get(h, i, a) / (2.0 * (moFb_->get(h, i, i) - moFb_->get(h, a, a)) + orbital_level_shift_);
                X_b_->set(h, i, a, value);
                X_b_->set(h, a, i, (-1.0) * value);
            }
        }
    }

    // Determine the rotation generator with respect to the reference orbitals
    Xtotal_a_->add(X_a_);
    Xtotal_b_->add(X_b_);

    dcft_timer_off("DCFTSolver::compute_orbital_rotation_jacobi()");
}

void
DCFTSolver::rotate_orbitals()
{
    dcft_timer_on("DCFTSolver::rotate_orbitals()");

    // Initialize the orbital rotation matrix
    SharedMatrix U_a(new Matrix("Orbital rotation matrix (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedMatrix U_b(new Matrix("Orbital rotation matrix (Beta)", nirrep_, nmopi_, nmopi_));

    // Compute the orbital rotation matrix and rotate the orbitals

    // U = I
    U_a->identity();
    U_b->identity();

    // U += X
    U_a->add(Xtotal_a_);
    U_b->add(Xtotal_b_);

    // U += 0.5 * X * X
    U_a->gemm(false, false, 0.5, Xtotal_a_, Xtotal_a_, 1.0);
    U_b->gemm(false, false, 0.5, Xtotal_b_, Xtotal_b_, 1.0);

    // Orthogonalize the U vectors
    int rowA = U_a->nrow();
    int colA = U_a->ncol();

    double **U_a_block = block_matrix(rowA, colA);
    memset(U_a_block[0], 0, sizeof(double)*rowA*colA);
    U_a_block = U_a->to_block_matrix();
    schmidt(U_a_block, rowA, colA, "outfile");
    U_a->set(U_a_block);
    free_block(U_a_block);

    int rowB = U_b->nrow();
    int colB = U_b->ncol();

    double **U_b_block = block_matrix(rowB, colB);
    memset(U_b_block[0], 0, sizeof(double)*rowB*colB);
    U_b_block = U_b->to_block_matrix();
    schmidt(U_b_block, rowB, colB, "outfile");
    U_b->set(U_b_block);
    free_block(U_b_block);

    // Rotate the orbitals
    Ca_->gemm(false, false, 1.0, old_ca_, U_a, 0.0);
    Cb_->gemm(false, false, 1.0, old_cb_, U_b, 0.0);

    dcft_timer_off("DCFTSolver::rotate_orbitals()");

}

}} //End namespaces
