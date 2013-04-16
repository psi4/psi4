#include <libtrans/integraltransform.h>
#include <libpsio/psio.hpp>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libdiis/diismanager.h>
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

void
DCFTSolver::run_simult_dcft_oo()
{
    // This is the simultaneous orbital/lambda update algorithm for the orbital-optimized methods
    int cycle = 0;

    fprintf(outfile, "\n\n\t*=================================================================================*\n"
                         "\t* Cycle   Max Orb Grad    RMS Lambda Error   delta E        Total Energy     DIIS *\n"
                         "\t*---------------------------------------------------------------------------------*\n");

    // Copy the reference orbitals and to use them as the reference for the orbital rotation
    old_ca_->copy(Ca_);
    old_cb_->copy(Cb_);

    // Set up the DIIS manager
    DIISManager diisManager(maxdiis_, "DCFT DIIS vectors");
    dpdbuf4 Laa, Lab, Lbb;
    dpd_buf4_init(&Laa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&Lbb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
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
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Lbb);

    while((!orbitalsDone_ || !cumulantDone_ || !densityConverged_ || !energyConverged_)
            && cycle++ < maxiter_){
        std::string diisString;
        // Save the old energy
        old_total_energy_ = new_total_energy_;
        // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
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
        new_total_energy_ = lambda_energy_;
        // Build new Tau from the density cumulant in the MO basis and transform it the SO basis
        build_tau();
        if (exact_tau_) {
            refine_tau();
        }
        transform_tau();
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
        // Compute new SCF energy
        compute_scf_energy();
        // Add SCF energy contribution to the total DCFT energy
        new_total_energy_ += scf_energy_;
        // Compute orbital gradient and check convergence
        orbitals_convergence_ = compute_orbital_residual();
        orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
        // Check convergence of the total DCFT energy
        energyConverged_ = fabs(old_total_energy_ - new_total_energy_) < cumulant_threshold_;
        // Compute the orbital rotation step using Jacobi method
        compute_orbital_rotation_jacobi();
        if(orbitals_convergence_ < diis_start_thresh_ && cumulant_convergence_ < diis_start_thresh_){
            //Store the DIIS vectors
            dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
            dpd_buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
            dpd_buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
            dpd_buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
            dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
            dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
            dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
            if(diisManager.add_entry(10, orbital_gradient_a_.get(), orbital_gradient_b_.get(), &Raa, &Rab, &Rbb,
                                     Xtotal_a_.get(), Xtotal_b_.get(), &Laa, &Lab, &Lbb)){
                diisString += "S";
            }
            if(diisManager.subspace_size() > mindiisvecs_){
                diisString += "/E";
                diisManager.extrapolate(5, Xtotal_a_.get(), Xtotal_b_.get(), &Laa, &Lab, &Lbb);
            }
            dpd_buf4_close(&Raa);
            dpd_buf4_close(&Rab);
            dpd_buf4_close(&Rbb);
            dpd_buf4_close(&Laa);
            dpd_buf4_close(&Lab);
            dpd_buf4_close(&Lbb);
        }
        // Obtain new orbitals
        rotate_orbitals();
        // Make sure that the orbital phase is retained
        if(!correct_mo_phases(false)){
            fprintf(outfile,"\t\tThere was a problem correcting the MO phases.\n"
                            "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
        }
        // Write orbitals to the checkpoint file
        write_orbitals_to_checkpoint();
        // Transform two-electron integrals to the MO basis using new orbitals, build denominators
        transform_integrals();
        // Update SCF density (Kappa) and check its RMS
        densityConverged_ = update_scf_density() < orbitals_threshold_;
        // If we've performed enough lambda updates since the last orbitals
        // update, reset the counter so another SCF update is performed
        fprintf(outfile, "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                cycle, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                new_total_energy_, diisString.c_str());
        fflush(outfile);
    }

    fprintf(outfile, "\t*=================================================================================*\n");
}

double
DCFTSolver::compute_orbital_residual() {

    dpdfile2 Xai, Xia;

    // Compute the unrelaxed densities for the orbital gradient
    compute_unrelaxed_density_OOOO();
    compute_unrelaxed_density_OOVV();
    compute_unrelaxed_density_OVOV();

    // Compute the OV part of the orbital gradient
    compute_orbital_gradient_OV();

    // Compute the VO part of the orbital gradient
    compute_orbital_gradient_VO();

    dpd_file2_init(&Xia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&Xai, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_mat_init(&Xia);
    dpd_file2_mat_init(&Xai);
    dpd_file2_mat_rd(&Xia);
    dpd_file2_mat_rd(&Xai);

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

    dpd_file2_close(&Xai);
    dpd_file2_close(&Xia);

    dpd_file2_init(&Xia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&Xai, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_file2_mat_init(&Xia);
    dpd_file2_mat_init(&Xai);
    dpd_file2_mat_rd(&Xia);
    dpd_file2_mat_rd(&Xai);

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

    dpd_file2_close(&Xai);
    dpd_file2_close(&Xia);

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
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_contract222(&H, &T, &X, 0, 1, 1.0, 0.0);
    dpd_file2_close(&T);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_ia = H_ib Tau_ba
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
    dpd_contract222(&H, &T, &X, 0, 1, 1.0, 0.0);
    dpd_file2_close(&T);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_OV: Two-electron contributions

    //
    // 2 * <OV||VV> Г_VVVV
    //

    // Compute contributions from VVVV density
    // 1. X_ia <-- <ib||cd> tau_ca tau_db
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    // Alpha contribution X_IA
    dpd_file2_init(&Y2_OV, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Y2 <O|V>");

    // Y2_IA = <IA|CD> Tau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    dpd_contract422(&I, &T_VV, &Y2_OV, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y2_IA -= (IA|CD) Tau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    dpd_contract422(&I, &T_VV, &Y2_OV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y2_IA -= (IA|cd) Tau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    dpd_contract422(&I, &T_vv, &Y2_OV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);

    // X_IA -= Y2_IC Tau_CA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_contract222(&Y2_OV, &T_VV, &X, 0, 1, -1.0, 1.0);
    dpd_file2_close(&X);

    dpd_file2_close(&Y2_OV);

    // Beta contribution X_ia
    dpd_file2_init(&Y2_ov, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Y2 <o|v>");

    // Y2_ib = <ib|cd> Tau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
    dpd_contract422(&I, &T_vv, &Y2_ov, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y2_ib -= (ib|cd) Tau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
    dpd_contract422(&I, &T_vv, &Y2_ov, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y2_ib -= (ib|CD) Tau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[V>=V]+"),
                  ID("[o,v]"), ID("[V>=V]+"), 0, "MO Ints (ov|VV)");
    dpd_contract422(&I, &T_VV, &Y2_ov, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);

    // X_ia -= Y2_ic Tau_ca
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_contract222(&Y2_ov, &T_vv, &X, 0, 1, -1.0, 1.0);
    dpd_file2_close(&X);

    dpd_file2_close(&Y2_ov);

    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    // 2. X_ia <-- 1/4 <ib||cd> lambda_abkl lambda_klcd

    //  W_IBKL = <IB||CD> lambda_KLCD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V>V]-"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    dpd_contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);


    //  W_KlIb = 2 lambda_KlCd <Ib|Cd>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    dpd_contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);

    //  W_LkBi = 2 lambda_LkDc <Bi|Dc>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    dpd_contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);

    //  W_ibkl = <ib||cd> lambda_klcd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v>v]-"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    dpd_contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);


    // X_IA +=  1/4 W_IBKL lambda_KLAB
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    dpd_buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");

    dpd_contract442(&W, &LL, &X, 0, 2, 0.25, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&LL);
    dpd_file2_close(&X);

    // X_IA +=  1/2 W_KlIb lambda_KlAb
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    dpd_buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

    dpd_contract442(&W, &LL, &X, 2, 2, 0.5, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&LL);
    dpd_file2_close(&X);

    // X_ia +=  1/2 W_LkBi lambda_LkBa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    dpd_buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

    dpd_contract442(&W, &LL, &X, 3, 3, 0.5, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&LL);
    dpd_file2_close(&X);

    // X_ia += 1/4 W_ibkl lambda_klab
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    dpd_buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");

    dpd_contract442(&W, &LL, &X, 0, 2, 0.25, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&LL);
    dpd_file2_close(&X);

    //
    // <OO||OV> Г_OOVV
    //

    // X_IA += <BI||JK> Г_BAJK
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V>V]-"), ID("[O>O]-"), 0, "Gamma <VV|OO>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_IA += 2 * <Ib|Jk> Г_AbJk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia += <bi||jk> Г_bajk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v>v]-"), ID("[o>o]-"), 0, "Gamma <vv|oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia += 2 * <Bi|Jk> Г_BaJk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <OO||OV> Г_OVOV
    //

    // X_IA += <JB||KI> Г_JBKA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 1, "MO Ints <OV|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_IA += <kB|jI> Г_kBjA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|oO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_IA -= <Kb|jI> Г_KbjA
    // Note: <Kb|jI> integrals are resorted <bK|jI> integrals.
    // <Kb||jI> Г_KbjA = (-1) * <bK|jI> * Г_KbjA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,O]"),
                  ID("[O,v]"), ID("[o,O]"), 0, "MO Ints <Ov|oO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia += <jb||ki> Г_jbka
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 1, "MO Ints <ov|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia += <Kb|Ji> Г_KbJa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia -= <kB|Ji> Г_kBJa
    // Note: <kB|Ji> integrals are resorted <Bk|Ji> integrals.
    // <kB||Ji> Г_kBJa = (-1) * <Bk|Ji> * Г_kBJa

    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[O,o]"),
                  ID("[o,V]"), ID("[O,o]"), 0, "MO Ints <oV|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "Gamma <oV|Ov>");

    dpd_contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

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
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    dpd_file2_mat_init(&X);
    dpd_file2_mat_init(&H);
    dpd_file2_mat_init(&T);
    dpd_file2_mat_rd(&H);
    dpd_file2_mat_rd(&T);
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
    dpd_file2_mat_wrt(&X);
    dpd_file2_close(&T);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_ai = H_ja Tau_ji
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    dpd_file2_mat_init(&X);
    dpd_file2_mat_init(&H);
    dpd_file2_mat_init(&T);
    dpd_file2_mat_rd(&H);
    dpd_file2_mat_rd(&T);
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
    dpd_file2_mat_wrt(&X);
    dpd_file2_close(&T);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_VO: Two-electron contributions

    //
    // 2 * <VO||OO> Г_OOOO
    //

    // X_AI += 2 * <AJ||KL> Г_IJKL
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>O]-"), ID("[O>O]-"), 0, "Gamma <OO|OO>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI += 4 * <Aj|Kl> Г_IjKl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += 2 * <aj||kl> Г_ijkl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>o]-"), ID("[o>o]-"), 0, "Gamma <oo|oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI += 4 * <Ja|Kl> Г_JiKl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <VO||VV> Г_OOVV
    //

    // X_AI += <JA||BC> Г_JIBC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Gamma <OO|VV>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI += <Aj|Bc> Г_IjBc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += <ja||bc> Г_jibc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Gamma <oo|vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += <Ja|Bc> Г_JiBc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <OV||VV> Г_OVOV
    //

    // X_AI += <JB||AC> Г_JBIC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI += <Jb|Ac> Г_JbIc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI -= <jB|Ac> Г_jBIc
    // Note: <jB|Ac> integrals are resorted <Bj|Ac> integrals.
    // <jB||Ac> Г_jBIc = (-1) * <Bj|Ac> Г_jBIc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[V,v]"),
                  ID("[o,V]"), ID("[V,v]"), 0, "MO Ints <oV|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "Gamma <oV|Ov>");

    dpd_contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += <jb||ac> Г_jbic
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += <jB|aC> Г_jBiC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai -= <Jb|aC> Г_JbiC
    // Note: <Jb|aC> integrals are resorted <bJ|aC> integrals.
    // <Jb||aC> Г_JbiC = (-1) * <bJ|aC> Г_JbiC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[v,V]"),
                  ID("[O,v]"), ID("[v,V]"), 0, "MO Ints <Ov|vV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void
DCFTSolver::compute_orbital_rotation_jacobi() {

    // Determine the orbital rotation step
    // Alpha spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = naoccpi_[h]; a < nmopi_[h]; ++a){
                double value = orbital_gradient_a_->get(h, i, a) / (2.0 * (moFa_->get(h, i, i) - moFa_->get(h, a, a)));
                X_a_->set(h, i, a, value);
                X_a_->set(h, a, i, (-1.0) * value);
            }
        }
    }

    // Beta spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = nboccpi_[h]; a < nmopi_[h]; ++a){
                double value = orbital_gradient_b_->get(h, i, a) / (2.0 * (moFb_->get(h, i, i) - moFb_->get(h, a, a)));
                X_b_->set(h, i, a, value);
                X_b_->set(h, a, i, (-1.0) * value);
            }
        }
    }

    // Determine the rotation generator with respect to the reference orbitals
    Xtotal_a_->add(X_a_);
    Xtotal_b_->add(X_b_);

}

void
DCFTSolver::rotate_orbitals()
{

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
    schmidt(U_a_block, rowA, colA, outfile);
    U_a->set(U_a_block);
    free_block(U_a_block);

    int rowB = U_b->nrow();
    int colB = U_b->ncol();

    double **U_b_block = block_matrix(rowB, colB);
    memset(U_b_block[0], 0, sizeof(double)*rowB*colB);
    U_b_block = U_b->to_block_matrix();
    schmidt(U_b_block, rowB, colB, outfile);
    U_b->set(U_b_block);
    free_block(U_b_block);

    // Rotate the orbitals
    Ca_->gemm(false, false, 1.0, old_ca_, U_a, 0.0);
    Cb_->gemm(false, false, 1.0, old_cb_, U_b, 0.0);


}

}} //End namespaces

