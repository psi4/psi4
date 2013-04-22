#include <libtrans/integraltransform.h>
#include <libpsio/psio.hpp>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libdiis/diismanager.h>
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

void
DCFTSolver::run_qc_dcft()
{
    // Quadratically-convergent algorithm: solution of the Newton-Raphson equations
    // for the simultaneous optimization of the cumulant and the orbitals

    if (options_.get_str("QC_TYPE") == "SIMULTANEOUS") {
        fprintf(outfile, "\n\n\t*==========================================================================================*\n"
                             "\t* Cycle   RMS Orb Grad   RMS Lambda Error    delta E         Total Energy     NI(NR)  DIIS *\n"
                             "\t*------------------------------------------------------------------------------------------*\n");
    }
    else {
        fprintf(outfile, "\n\n\t*====================================================================================================*\n"
                             "\t* Cycle   RMS Orb Grad   RMS Lambda Error    delta E         Total Energy     NI(Orb)  NI(Cum)  DIIS *\n"
                             "\t*----------------------------------------------------------------------------------------------------*\n");
    }

    int cycle = 0;
    int cycle_NR = 0;
    int cycle_jacobi = 0;

    // Copy the reference orbitals and to use them as the reference for the orbital rotation
    old_ca_->copy(Ca_);
    old_cb_->copy(Cb_);

    orbitals_convergence_ = compute_scf_error_vector();

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

    while((!orbitalsDone_ || !cumulantDone_ || !energyConverged_ || !densityConverged_) && cycle++ < maxiter_ ) {

        std::string diisString;
        // Compute the generalized Fock matrix and orbital gradient in the MO basis
        compute_orbital_gradient();
        // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
        build_cumulant_intermediates();
        // Compute the residuals for density cumulant equations
        cumulant_convergence_ = compute_cumulant_residual();
        // Save the old energy
        old_total_energy_ = new_total_energy_;
        // Compute new SCF energy
        compute_scf_energy();
        // Add SCF energy contribution to the total DCFT energy
        new_total_energy_ = scf_energy_;
        // Compute new DCFT energy (lambda contribution)
        compute_dcft_energy();
        // Add lambda energy to the DCFT total energy
        new_total_energy_ += lambda_energy_;
        // Check convergence of the total DCFT energy
        energyConverged_ = fabs(old_total_energy_ - new_total_energy_) < cumulant_threshold_;
        // Determine the independent pairs (IDPs) and create array for the orbital and cumulant gradient in the basis of IDPs
        form_idps();
        if (nidp_ != 0) {
            // Compute sigma vector in the basis of IDPs
            compute_sigma_vector();
            // Solve the NR equations using conjugate gradients
            cycle_NR = iterate_nr_conjugate_gradients();
            // Check the convergence by computing the change in the orbitals and the cumulant
            check_qc_convergence();
            orbitalsDone_ = orbitals_convergence_ < orbitals_threshold_;
            cumulantDone_ = cumulant_convergence_ < cumulant_threshold_;
            // Update cumulant first
            if (options_.get_str("QC_TYPE") == "SIMULTANEOUS") {
                update_cumulant_nr();
            }
            else {
                double old_energy = old_total_energy_;
                cycle_jacobi = run_twostep_dcft_cumulant_updates();
                old_total_energy_ = old_energy;
            }
            // Compute the rotation for the orbitals
            compute_orbital_rotation_nr();
            // DIIS
            if(orbitals_convergence_ < diis_start_thresh_ && cumulant_convergence_ < diis_start_thresh_){
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
            // Update orbitals
            rotate_orbitals();
            // Print the iterative trace
            if (options_.get_str("QC_TYPE") == "SIMULTANEOUS") {
                fprintf(outfile, "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f   %3d    %-3s *\n",
                        cycle, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                        new_total_energy_, cycle_NR, diisString.c_str());
            }
            else {
                fprintf(outfile, "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f   %3d      %3d      %-3s *\n",
                        cycle, orbitals_convergence_, cumulant_convergence_, new_total_energy_ - old_total_energy_,
                        new_total_energy_, cycle_NR, cycle_jacobi, diisString.c_str());

            }
            fflush(outfile);
            if (orbital_idp_ != 0) {
                // Update the density
                densityConverged_ = update_scf_density() < orbitals_threshold_;
                // Write orbitals to the checkpoint file
                write_orbitals_to_checkpoint();
                // Transform two-electron integrals to the MO basis using new orbitals, build denominators
                // TODO: Transform_integrals shouldn't call build denominators for the QC alogorithm
                transform_integrals();
            }
        }
        else break;

    }

    if (options_.get_str("QC_TYPE") == "SIMULTANEOUS") {
        fprintf(outfile, "\t*==========================================================================================*\n");
    }
    else {
        fprintf(outfile, "\t*====================================================================================================*\n");
    }

    if(!orbitalsDone_ || !cumulantDone_ || !densityConverged_ || !energyConverged_)
        throw ConvergenceError<int>("DCFT", maxiter_, cumulant_threshold_,
                               cumulant_convergence_, __FILE__, __LINE__);

}

void
DCFTSolver::compute_orbital_gradient(){

    // Build guess Tau from the density cumulant in the MO basis and transform it to the SO basis
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
    // Form F0 matrix
    moF0a_->copy(Fa_);
    moF0b_->copy(Fb_);
    // Transform the F0 matrix to the MO basis
    moF0a_->transform(Ca_);
    moF0b_->transform(Cb_);
    // Add non-idempotent density contribution (Tau) to the Fock matrix: F += Gbar * Tau
    Fa_->add(g_tau_a_);
    Fb_->add(g_tau_b_);
    // Copy the SO basis Fock for the transformation to the MO basis
    moFa_->copy(Fa_);
    moFb_->copy(Fb_);
    // Transform the Fock matrix to the MO basis
    moFa_->transform(Ca_);
    moFb_->transform(Cb_);
    // Update the Fock matrix DPD file2
    dpdfile2 F;
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    //Alpha Occupied
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_file2_mat_init(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j < naoccpi_[h]; ++j){
                F.matrix[h][i][j] = moFa_->get(h, i, j);
            }
        }
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    //Alpha Virtual
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    dpd_file2_mat_init(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < navirpi_[h]; ++i){
            for(int j = 0 ; j < navirpi_[h]; ++j){
                F.matrix[h][i][j] = moFa_->get(h, i + naoccpi_[h], j + naoccpi_[h]);
            }
        }
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    //Beta Occupied
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_file2_mat_init(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j < nboccpi_[h]; ++j){
                F.matrix[h][i][j] = moFb_->get(h, i, j);
            }
        }
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    //Beta Virtual
    dpd_file2_init(&F, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    dpd_file2_mat_init(&F);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nbvirpi_[h]; ++i){
            for(int j = 0 ; j < nbvirpi_[h]; ++j){
                F.matrix[h][i][j] = moFb_->get(h, i + nboccpi_[h], j + nboccpi_[h]);
            }
        }
    }
    dpd_file2_mat_wrt(&F);
    dpd_file2_close(&F);

    psio_->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // Initialize the idempotent contribution to the OPDM (Kappa)
    if (!orbital_optimized_) {
        SharedMatrix full_kappa_a(new Matrix("MO basis Full Kappa (Alpha)", nirrep_, nmopi_, nmopi_));
        SharedMatrix full_kappa_b(new Matrix("MO basis Full Kappa (Beta)", nirrep_, nmopi_, nmopi_));
        // Compute Kappa in the MO basis
        for(int h = 0; h < nirrep_; ++h){
            for(int i = 0; i < naoccpi_[h]; ++i){
                full_kappa_a->set(h, i, i, 1.0);
            }
            for(int i = 0; i < nboccpi_[h]; ++i){
                full_kappa_b->set(h, i, i, 1.0);
            }
        }

        // Form F * kappa (Alpha spin)
        orbital_gradient_a_->gemm(false, false, 2.0, moFa_, full_kappa_a, 0.0);
        // Form -kappa * F (Alpha spin) to obtain [F,kappa]
        orbital_gradient_a_->gemm(false, false, -2.0, full_kappa_a, moFa_, 1.0);

        // Form F * kappa (Beta spin)
        orbital_gradient_b_->gemm(false, false, 2.0, moFb_, full_kappa_b, 0.0);
        // Form -kappa * F (Alpha spin) to obtain [F,kappa]
        orbital_gradient_b_->gemm(false, false, -2.0, full_kappa_b, moFb_, 1.0);
    }
    else {
        compute_orbital_residual();
    }

}

void DCFTSolver::form_idps(){

    // Ignore orbital and cumulant gradient elements that are less than this value
    double cutoff = ((1.0e-10 < (cumulant_threshold_  * 0.01)) ? 1.0e-10 : (cumulant_threshold_  * 0.01));

    // Zero out the counters
    int old_nidp = nidp_;
    nidp_ = 0;
    orbital_idp_a_ = 0;
    orbital_idp_b_ = 0;
    orbital_idp_ = 0;
    cumulant_idp_aa_ = 0;
    cumulant_idp_ab_ = 0;
    cumulant_idp_bb_ = 0;
    cumulant_idp_ = 0;

    // Zero the lookup array
    ::memset(lookup_orbitals_, '\0', sizeof(int)*dim_orbitals_);

    // Temporary vectors containing gradient value and diagonal part of the Hessian for each IDP
    double *grad = new double[dim_];
    double *Hd = new double[dim_];

    // Count the number of IDPs for orbital rotations (Alpha spin)
    // The minus sign in the gradient takes into account the sign of the g vector in the N-R equations: dX H = -g
    int orbital_address = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                if (fabs(orbital_gradient_a_->get(h, i, a + naoccpi_[h])) > cutoff) {
                    lookup_orbitals_[orbital_address] = 1;
                    grad[orbital_idp_a_] = (-1.0) * orbital_gradient_a_->get(h, i, a + naoccpi_[h]);
                    Hd[orbital_idp_a_] = 2.0 * (moFa_->get(h, a + naoccpi_[h], a + naoccpi_[h]) - moFa_->get(h, i, i));
                    orbital_idp_a_++;
                }
                orbital_address++;
            }
        }
    }

    // Count the number of IDPs for orbital rotations (Beta spin)
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = 0; a < nbvirpi_[h]; ++a){
                if (fabs(orbital_gradient_b_->get(h, i, a + nboccpi_[h])) > cutoff) {
                    lookup_orbitals_[orbital_address] = 1;
                    int index = orbital_idp_a_ + orbital_idp_b_;
                    grad[index] = (-1.0) * orbital_gradient_b_->get(h, i, a + nboccpi_[h]);
                    Hd[index] = 2.0 * (moFb_->get(h, a + nboccpi_[h], a + nboccpi_[h]) - moFb_->get(h, i, i));
                    orbital_idp_b_++;
                }
                orbital_address++;
            }
        }
    }

    orbital_idp_ = orbital_idp_a_ + orbital_idp_b_;

    if(options_.get_str("QC_TYPE") == "SIMULTANEOUS") {
        // Count the number of IDPs for cumulant updates
        dpdbuf4 R;

        ::memset(lookup_cumulant_, '\0', sizeof(int)*dim_cumulant_);

        int cumulant_address = 0;

        // Alpha-Alpha spin
        dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
        for(int h = 0; h < nirrep_; ++h){
            dpd_buf4_mat_irrep_init(&R, h);
            dpd_buf4_mat_irrep_rd(&R, h);
            #pragma omp parallel for
            for(size_t ij = 0; ij < R.params->rowtot[h]; ++ij){
                size_t i = R.params->roworb[h][ij][0];
                int isym = R.params->psym[i];
                i -= R.params->poff[isym];
                size_t j = R.params->roworb[h][ij][1];
                int jsym = R.params->qsym[j];
                j -= R.params->qoff[jsym];
                for(size_t ab = 0; ab < R.params->coltot[h]; ++ab){
                    size_t a = R.params->colorb[h][ab][0];
                    int asym = R.params->rsym[a];
                    a -= R.params->roff[asym];
                    size_t b = R.params->colorb[h][ab][1];
                    int bsym = R.params->ssym[b];
                    b -= R.params->soff[bsym];
                    if (fabs(R.matrix[h][ij][ab]) > cutoff) {
                        lookup_cumulant_[cumulant_address] = 1;
                        int index = orbital_idp_ + cumulant_idp_aa_;
                        grad[index] = -0.25 * R.matrix[h][ij][ab];
                        double value = moFa_->get(asym, a + naoccpi_[asym], a + naoccpi_[asym])
                                + moFa_->get(bsym, b + naoccpi_[bsym], b + naoccpi_[bsym])
                                - moFa_->get(isym, i, i) - moFa_->get(jsym, j, j);
                        Hd[index] = (1.0/16.0) * value;
                        cumulant_idp_aa_++ ;
                    }
                    cumulant_address++;
                }
            }
            dpd_buf4_mat_irrep_close(&R, h);
        }
        dpd_buf4_close(&R);

        // Alpha-Beta spin
        dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
        for(int h = 0; h < nirrep_; ++h){
            dpd_buf4_mat_irrep_init(&R, h);
            dpd_buf4_mat_irrep_rd(&R, h);
            #pragma omp parallel for
            for(size_t ij = 0; ij < R.params->rowtot[h]; ++ij){
                size_t i = R.params->roworb[h][ij][0];
                int isym = R.params->psym[i];
                i -= R.params->poff[isym];
                size_t j = R.params->roworb[h][ij][1];
                int jsym = R.params->qsym[j];
                j -= R.params->qoff[jsym];
                for(size_t ab = 0; ab < R.params->coltot[h]; ++ab){
                    size_t a = R.params->colorb[h][ab][0];
                    int asym = R.params->rsym[a];
                    a -= R.params->roff[asym];
                    size_t b = R.params->colorb[h][ab][1];
                    int bsym = R.params->ssym[b];
                    b -= R.params->soff[bsym];
                    if (fabs(R.matrix[h][ij][ab]) > cutoff) {
                        lookup_cumulant_[cumulant_address] = 1;
                        int index = orbital_idp_ + cumulant_idp_aa_ + cumulant_idp_ab_;
                        grad[index] = -0.25 * R.matrix[h][ij][ab];
                        double value = moFa_->get(asym, a + naoccpi_[asym], a + naoccpi_[asym])
                                + moFb_->get(bsym, b + nboccpi_[bsym], b + nboccpi_[bsym])
                                - moFa_->get(isym, i, i) - moFb_->get(jsym, j, j);
                        Hd[index] = (1.0/16.0) * value;
                        cumulant_idp_ab_++ ;
                    }
                    cumulant_address++;
                }
            }
            dpd_buf4_mat_irrep_close(&R, h);
        }
        dpd_buf4_close(&R);

        // Beta-Beta spin
        dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
        for(int h = 0; h < nirrep_; ++h){
            dpd_buf4_mat_irrep_init(&R, h);
            dpd_buf4_mat_irrep_rd(&R, h);
            #pragma omp parallel for
            for(size_t ij = 0; ij < R.params->rowtot[h]; ++ij){
                size_t i = R.params->roworb[h][ij][0];
                int isym = R.params->psym[i];
                i -= R.params->poff[isym];
                size_t j = R.params->roworb[h][ij][1];
                int jsym = R.params->qsym[j];
                j -= R.params->qoff[jsym];
                for(size_t ab = 0; ab < R.params->coltot[h]; ++ab){
                    size_t a = R.params->colorb[h][ab][0];
                    int asym = R.params->rsym[a];
                    a -= R.params->roff[asym];
                    size_t b = R.params->colorb[h][ab][1];
                    int bsym = R.params->ssym[b];
                    b -= R.params->soff[bsym];
                    if (fabs(R.matrix[h][ij][ab]) > cutoff) {
                        lookup_cumulant_[cumulant_address] = 1;
                        int index = orbital_idp_ + cumulant_idp_aa_ + cumulant_idp_ab_ + cumulant_idp_bb_;
                        grad[index] = -0.25 * R.matrix[h][ij][ab];
                        double value = moFb_->get(asym, a + nboccpi_[asym], a + nboccpi_[asym])
                                + moFb_->get(bsym, b + nboccpi_[bsym], b + nboccpi_[bsym])
                                - moFb_->get(isym, i, i) - moFb_->get(jsym, j, j);
                        Hd[index] = (1.0/16.0) * value;
                        cumulant_idp_bb_++ ;
                    }
                    cumulant_address++;
                }
            }
            dpd_buf4_mat_irrep_close(&R, h);
        }
        dpd_buf4_close(&R);

        cumulant_idp_ = cumulant_idp_aa_ + cumulant_idp_ab_ + cumulant_idp_bb_;
    }

    nidp_ = orbital_idp_ + cumulant_idp_;

    // Reallocate the memory for the arrays if the dimensions changed
    if (old_nidp != nidp_) {
        gradient_ = SharedVector(new Vector("Orbital and Cumulant Gradient in the IDP basis", nidp_));
        Hd_ = SharedVector(new Vector("Diagonal part of the Hessian in the IDP basis", nidp_));
        X_ = SharedVector(new Vector("Step vector in the IDP basis", nidp_));
        sigma_ = SharedVector(new Vector("Sigma vector in the IDP basis", nidp_));
        D_ = SharedVector(new Vector("Conjugate direction vector in the IDP basis", nidp_));
        R_ = SharedVector(new Vector("Residual vector in the IDP basis", nidp_));
        S_ = SharedVector(new Vector("Search direction vector in the IDP basis", nidp_));
        Q_ = SharedVector(new Vector("New element of the Krylov subspace vector in the IDP basis", nidp_));
    }

    R_->zero();
    S_->zero();
    Q_->zero();

    // Set the gradient Hessian diagonal vectors
    for (int p = 0; p < nidp_; ++p) {
        gradient_->set(p, grad[p]);
        Hd_->set(p, Hd[p]);
        // Compute the guess for the step vector
        X_->set(p, grad[p]/Hd[p]);
    }

    // Perform stability analysis if requested by the user
    if (options_.get_bool("STABILITY_CHECK")) run_davidson();

    sigma_->zero();

    // Set the gradient Hessian diagonal vectors
    for (int p = 0; p < nidp_; ++p) {
        // D is used to store X for the initial formation of the sigma vector
        D_->set(p, X_->get(p));
    }

    delete[] Hd;
    delete[] grad;

}

void
DCFTSolver::compute_sigma_vector() {

    dpdfile2 D2, S2;
    dpdbuf4 D4, S4;

    // Copy the IDP arrays from memory to DPD file2 and buf4
    int orbital_address = 0;
    int idpcount = 0;

    // Conjugate direction for orbital rotations (Alpha spin)
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_file2_mat_init(&D2);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                D2.matrix[h][i][a] = 0.0;
                if (lookup_orbitals_[orbital_address]) {
                    D2.matrix[h][i][a] = D_->get(idpcount);
                    idpcount++;
                }
                orbital_address++;
            }
        }
    }
    dpd_file2_mat_wrt(&D2);
    dpd_file2_close(&D2);

    // Conjugate direction for orbital rotations (Beta spin)
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_file2_mat_init(&D2);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = 0; a < nbvirpi_[h]; ++a){
                D2.matrix[h][i][a] = 0.0;
                if (lookup_orbitals_[orbital_address]) {
                    D2.matrix[h][i][a] = D_->get(idpcount);
                    idpcount++;
                }
                orbital_address++;
            }
        }
    }
    dpd_file2_mat_wrt(&D2);
    dpd_file2_close(&D2);

    if(options_.get_str("QC_TYPE") == "SIMULTANEOUS") {

        int cumulant_address = 0;

        // Conjugate directions for cumulant updates
        // Alpha-Alpha spin
        dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "D4 <OO|VV>");
        for(int h = 0; h < nirrep_; ++h){
            dpd_buf4_mat_irrep_init(&D4, h);
            #pragma omp parallel for
            for(size_t ij = 0; ij < D4.params->rowtot[h]; ++ij){
                for(size_t ab = 0; ab < D4.params->coltot[h]; ++ab){
                    D4.matrix[h][ij][ab] = 0.0;
                    if (lookup_cumulant_[cumulant_address]) {
                        D4.matrix[h][ij][ab] = D_->get(idpcount);
                        idpcount++;
                    }
                    cumulant_address++;
                }
            }
            dpd_buf4_mat_irrep_wrt(&D4, h);
            dpd_buf4_mat_irrep_close(&D4, h);
        }
        dpd_buf4_close(&D4);

        // Alpha-Beta spin
        dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "D4 <Oo|Vv>");
        for(int h = 0; h < nirrep_; ++h){
            dpd_buf4_mat_irrep_init(&D4, h);
            #pragma omp parallel for
            for(size_t ij = 0; ij < D4.params->rowtot[h]; ++ij){
                for(size_t ab = 0; ab < D4.params->coltot[h]; ++ab){
                    D4.matrix[h][ij][ab] = 0.0;
                    if (lookup_cumulant_[cumulant_address]) {
                        D4.matrix[h][ij][ab] = D_->get(idpcount);
                        idpcount++;
                    }
                    cumulant_address++;
                }
            }
            dpd_buf4_mat_irrep_wrt(&D4, h);
            dpd_buf4_mat_irrep_close(&D4, h);
        }
        dpd_buf4_close(&D4);

        // Beta-Beta spin
        dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "D4 <oo|vv>");
        for(int h = 0; h < nirrep_; ++h){
            dpd_buf4_mat_irrep_init(&D4, h);
            #pragma omp parallel for
            for(size_t ij = 0; ij < D4.params->rowtot[h]; ++ij){
                for(size_t ab = 0; ab < D4.params->coltot[h]; ++ab){
                    D4.matrix[h][ij][ab] = 0.0;
                    if (lookup_cumulant_[cumulant_address]) {
                        D4.matrix[h][ij][ab] = D_->get(idpcount);
                        idpcount++;
                    }
                    cumulant_address++;
                }
            }
            dpd_buf4_mat_irrep_wrt(&D4, h);
            dpd_buf4_mat_irrep_close(&D4, h);
        }
        dpd_buf4_close(&D4);
    }

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // Orbital-Orbital part of the Hessian
    compute_sigma_vector_orb_orb();
    // Orbital-Cumulant part of the Hessian
    if(options_.get_str("QC_TYPE") == "SIMULTANEOUS") {
        if (options_.get_bool("QC_COUPLING")) {
            compute_sigma_vector_orb_cum();
        }
        // Cumulant-Cumulant part of the Hessian
        compute_sigma_vector_cum_cum();
        // Cumulant-Orbital part of the Hessian
        if (options_.get_bool("QC_COUPLING")) {
            compute_sigma_vector_cum_orb();
        }
    }

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    // Now copy the sigma vector back to memory
    orbital_address = 0;
    idpcount = 0;
    // Alpha spin
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_mat_init(&S2);
    dpd_file2_mat_rd(&S2);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                if (lookup_orbitals_[orbital_address]) {
                    sigma_->set(idpcount, S2.matrix[h][i][a]);
                    idpcount++;
                }
                orbital_address++;
            }
        }
    }
    dpd_file2_close(&S2);

    // Beta spin
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_mat_init(&S2);
    dpd_file2_mat_rd(&S2);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = 0; a < nbvirpi_[h]; ++a){
                if (lookup_orbitals_[orbital_address]) {
                    sigma_->set(idpcount, S2.matrix[h][i][a]);
                    idpcount++;
                }
                orbital_address++;
            }
        }
    }
    dpd_file2_close(&S2);

    if(options_.get_str("QC_TYPE") == "SIMULTANEOUS") {

        int cumulant_address = 0;

        // Alpha-Alpha spin
        dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Sigma <OO|VV>");
        for(int h = 0; h < nirrep_; ++h){
            dpd_buf4_mat_irrep_init(&S4, h);
            dpd_buf4_mat_irrep_rd(&S4, h);
            #pragma omp parallel for
            for(size_t ij = 0; ij < S4.params->rowtot[h]; ++ij){
                for(size_t ab = 0; ab < S4.params->coltot[h]; ++ab){
                    if (lookup_cumulant_[cumulant_address]) {
                        sigma_->set(idpcount, S4.matrix[h][ij][ab]);
                        idpcount++;
                    }
                    cumulant_address++;
                }
            }
            dpd_buf4_mat_irrep_close(&S4, h);
        }
        dpd_buf4_close(&S4);

        // Alpha-Beta spin
        dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Sigma <Oo|Vv>");
        for(int h = 0; h < nirrep_; ++h){
            dpd_buf4_mat_irrep_init(&S4, h);
            dpd_buf4_mat_irrep_rd(&S4, h);
            #pragma omp parallel for
            for(size_t ij = 0; ij < S4.params->rowtot[h]; ++ij){
                for(size_t ab = 0; ab < S4.params->coltot[h]; ++ab){
                    if (lookup_cumulant_[cumulant_address]) {
                        sigma_->set(idpcount, S4.matrix[h][ij][ab]);
                        idpcount++;
                    }
                    cumulant_address++;
                }
            }
            dpd_buf4_mat_irrep_close(&S4, h);
        }
        dpd_buf4_close(&S4);

        // Beta-Beta spin
        dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Sigma <oo|vv>");
        for(int h = 0; h < nirrep_; ++h){
            dpd_buf4_mat_irrep_init(&S4, h);
            dpd_buf4_mat_irrep_rd(&S4, h);
            #pragma omp parallel for
            for(size_t ij = 0; ij < S4.params->rowtot[h]; ++ij){
                for(size_t ab = 0; ab < S4.params->coltot[h]; ++ab){
                    if (lookup_cumulant_[cumulant_address]) {
                        sigma_->set(idpcount, S4.matrix[h][ij][ab]);
                        idpcount++;
                    }
                    cumulant_address++;
                }
            }
            dpd_buf4_mat_irrep_close(&S4, h);
        }
        dpd_buf4_close(&S4);
    }

}

void
DCFTSolver::compute_sigma_vector_orb_orb() {

    dpdfile2 D2, S2;
    dpdbuf4 I;

    // Orbital-Orbital block of the Hessian

    // Sigma_IA = 2.0 * (2(IA|JB) - <AI|JB> - <IA|JB>) D_JB
    // Sigma_IA = 4.0 * (IA|JB) * D_JB
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    dpd_contract422(&I, &D2, &S2, 0, 0, 4.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&S2);

    // Sigma_IA += 4.0 * (IA|jb) * D_jb
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    dpd_contract422(&I, &D2, &S2, 0, 0, 4.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&S2);

    // Sigma_IA -= 2.0 * <IA|JB> * D_JB
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_contract422(&I, &D2, &S2, 0, 0, -2.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&S2);

    // Sigma_IA -= 2.0 * <AI|JB> * D_JB
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                  ID("[V,O]"), ID("[O,V]"), 0, "MO Ints <VO|OV>");
    dpd_contract422(&I, &D2, &S2, 0, 1, -2.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&S2);

    // Sigma_ia = 2.0 * (2(ia|jb) - <ai|jb> - <ia|jb>) D_jb
    // Sigma_ia = 4.0 * (ia|jb) * D_jb
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    dpd_contract422(&I, &D2, &S2, 0, 0, 4.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&S2);

    // Sigma_ia += 4.0 * (ia|JB) * D_JB
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[O,V]"),
                  ID("[o,v]"), ID("[O,V]"), 0, "MO Ints (ov|OV)");
    dpd_contract422(&I, &D2, &S2, 0, 0, 4.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&S2);

    // Sigma_ia -= 2.0 * <ia|jb> * D_jb
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
    dpd_contract422(&I, &D2, &S2, 0, 0, -2.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&S2);

    // Sigma_ia -= 2.0 * <ai|jb> * D_jb
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                  ID("[v,o]"), ID("[o,v]"), 0, "MO Ints <vo|ov>");
    dpd_contract422(&I, &D2, &S2, 0, 1, -2.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&S2);

}

void
DCFTSolver::compute_sigma_vector_orb_cum() {

    dpdfile2 S2, LD_OO, LD_oo, LD_VV, LD_vv;
    dpdbuf4 I, L, D4;

    // Orbital-Lambda block of the Hessian

    // Sigma_IA += 0.25 * (<AM||IK> + <AK||IM>) * Lambda_MLCD * D_KLCD
    // Sigma_IA -= 0.25 * (<IC||AE> + <IE||AC>) * Lambda_KLED * D_KLCD

    // Compute Lambda * D intermediates for both terms

    dpd_file2_init(&LD_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    dpd_file2_init(&LD_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");
    dpd_file2_init(&LD_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    dpd_file2_init(&LD_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");

    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "D4 <OO|VV>");

    // LD_IJ = Lambda_MLCD * D_KLCD
    dpd_contract442(&L, &D4, &LD_OO, 0, 0, 1.0, 0.0);

    // LD_AB = Lambda_KLED * D_KLCD
    dpd_contract442(&L, &D4, &LD_VV, 2, 2, 1.0, 0.0);

    dpd_buf4_close(&D4);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "D4 <oo|vv>");

    // LD_ij = Lambda_mlcd * D_klcd
    dpd_contract442(&L, &D4, &LD_oo, 0, 0, 1.0, 0.0);

    // LD_ab = Lambda_kled * D_klcd
    dpd_contract442(&L, &D4, &LD_vv, 2, 2, 1.0, 0.0);

    dpd_buf4_close(&D4);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D4 <Oo|Vv>");

    // LD_IJ += Lambda_IkAb D4_JkAb + Lambda_IkaB D4_JkaB
    dpd_contract442(&L, &D4, &LD_OO, 0, 0, 2.0, 1.0);
    // LD_ij += Lambda_KiAb D4_KjAb + Lambda_KiaB D4_KjaB
    dpd_contract442(&L, &D4, &LD_oo, 1, 1, 2.0, 1.0);
    // LD_AB += Lambda_IjAc D4_IjBc + Lambda_iJAc D4_iJBc
    dpd_contract442(&L, &D4, &LD_VV, 2, 2, 2.0, 1.0);
    // LD_ab += Lambda_IjCa D4_IjCb + Lambda_iJCa D4_iJCb
    dpd_contract442(&L, &D4, &LD_vv, 3, 3, 2.0, 1.0);

    dpd_buf4_close(&D4);
    dpd_buf4_close(&L);

    dpd_file2_close(&LD_OO);
    dpd_file2_close(&LD_oo);
    dpd_file2_close(&LD_VV);
    dpd_file2_close(&LD_vv);

    // Sigma_IA += 0.25 * (<AM||IK> + <AK||IM>) * LD_MK
    // Sigma_IA += 0.25 * (2(MK|AI) - <MK|IA> - <MK|AI>) LD_MK
    // Sigma_IA += 0.5 * (AI|MK) * LD_MK
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&LD_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O>=O]+"), 0, "MO Ints (VO|OO)");
    dpd_contract422(&I, &LD_OO, &S2, 0, 1, 0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_OO);
    dpd_file2_close(&S2);

    // Sigma_IA += 0.5 * (AI|mk) * LD_mk
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&LD_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,o]"),
                  ID("[V,O]"), ID("[o>=o]+"), 0, "MO Ints (VO|oo)");
    dpd_contract422(&I, &LD_oo, &S2, 0, 1, 0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_oo);
    dpd_file2_close(&S2);

    // Sigma_IA -= 0.25 * <IA|MK>  * LD_MK
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&LD_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "MO Ints <OV|OO>");
    dpd_contract422(&I, &LD_OO, &S2, 0, 0, -0.25, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_OO);
    dpd_file2_close(&S2);

    // Sigma_IA -= 0.25 * <AI|MK>  * LD_MK
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&LD_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <VO|OO>");
    dpd_contract422(&I, &LD_OO, &S2, 0, 1, -0.25, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_OO);
    dpd_file2_close(&S2);

    // Sigma_ia += 0.25 * (<am||ik> + <ak||im>) * LD_mk
    // Sigma_ia += 0.25 * (2(mk|ai) - <mk|ia> - <mk|ai>) LD_mk
    // Sigma_ia += 0.5 * (ai|mk) * LD_mk
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&LD_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o>=o]+"), 0, "MO Ints (vo|oo)");
    dpd_contract422(&I, &LD_oo, &S2, 0, 1, 0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_oo);
    dpd_file2_close(&S2);

    // Sigma_ia += 0.5 * (ai|MK) * LD_MK
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&LD_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[O,O]"),
                  ID("[v,o]"), ID("[O>=O]+"), 0, "MO Ints (vo|OO)");
    dpd_contract422(&I, &LD_OO, &S2, 0, 1, 0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_OO);
    dpd_file2_close(&S2);

    // Sigma_ia -= 0.25 * <ia|mk>  * LD_mk
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&LD_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 0, "MO Ints <ov|oo>");
    dpd_contract422(&I, &LD_oo, &S2, 0, 0, -0.25, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_oo);
    dpd_file2_close(&S2);

    // Sigma_ia -= 0.25 * <ai|mk>  * LD_mk
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&LD_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <vo|oo>");
    dpd_contract422(&I, &LD_oo, &S2, 0, 1, -0.25, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_oo);
    dpd_file2_close(&S2);

    // Sigma_IA -= 0.25 * (<IC||AE> + <IE||AC>) * Lambda_KLED * D_KLCD
    // Sigma_IA -= 0.25 * (2(EC|AI) - <IA|EC> - <AI|EC>) LD_EC
    // Sigma_IA -= 0.5 * (AI|CE) LD_CE
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&LD_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    dpd_contract422(&I, &LD_VV, &S2, 0, 0, -0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_VV);
    dpd_file2_close(&S2);

    // Sigma_IA -= 0.5 * (IA|ce) * LD_ce
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&LD_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    dpd_contract422(&I, &LD_vv, &S2, 0, 0, -0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_vv);
    dpd_file2_close(&S2);

    // Sigma_IA += 0.25 * <IA|CE> * LD_CE
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&LD_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    dpd_contract422(&I, &LD_VV, &S2, 1, 0, 0.25, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_VV);
    dpd_file2_close(&S2);

    // Sigma_IA += 0.25 * <IA|EC> * LD_EC
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_init(&LD_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    dpd_contract422(&I, &LD_VV, &S2, 0, 0, 0.25, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_VV);
    dpd_file2_close(&S2);

    // Sigma_ia -= 0.25 * (<ic||ae> + <ie||ac>) * Lambda_kled * D_klcd
    // Sigma_ia -= 0.25 * (2(ec|ai) - <ia|ec> - <ai|ec>) LD_ec
    // Sigma_ia -= 0.5 * (ia|ce) LD_ce
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&LD_VV, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
    dpd_contract422(&I, &LD_VV, &S2, 0, 0, -0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_VV);
    dpd_file2_close(&S2);

    // Sigma_ia -= 0.5 * (ia|CE) * LD_CE
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&LD_vv, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[V,V]"),
                  ID("[o,v]"), ID("[V>=V]+"), 0, "MO Ints (ov|VV)");
    dpd_contract422(&I, &LD_vv, &S2, 0, 0, -0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_vv);
    dpd_file2_close(&S2);

    // Sigma_ia += 0.25 * <ia|ce> * LD_ce
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&LD_VV, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
    dpd_contract422(&I, &LD_VV, &S2, 1, 0, 0.25, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_VV);
    dpd_file2_close(&S2);

    // Sigma_ia += 0.25 * <ia|ec> * LD_ec
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Sigma <o|v>");
    dpd_file2_init(&LD_VV, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
    dpd_contract422(&I, &LD_VV, &S2, 0, 0, 0.25, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&LD_VV);
    dpd_file2_close(&S2);

}

void
DCFTSolver::compute_sigma_vector_cum_cum() {

    dpdbuf4 I, D4, S4, T, Taa, Tab, Tbb, D4aa, D4ab, D4bb;

    // Lambda-Lambda block of the Hessian

    /*
     * Sigma_ijab += 1/16 Sum_cd gbar_cdab D_ijcd
     */
    // S_IJAB += (C>D) 1/16 Sum_CD gbar_CDAB D_IJCD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V>V]-"), ID("[V>V]-"),
                  ID("[V,V]"), ID("[V,V]"), 1, "MO Ints <VV|VV>");
    dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "D4 <OO|VV>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Sigma <OO|VV>");
    dpd_contract444(&D4, &I, &S4, 0, 0, 1.0/16.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&D4);
    dpd_buf4_close(&S4);

    // S_IjAb += 1/16 Sum_Cd g_CdAb D_IjCd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D4 <Oo|Vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Sigma <Oo|Vv>");
    dpd_contract444(&D4, &I, &S4, 0, 0, 1.0/16.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&D4);
    dpd_buf4_close(&S4);

    // S_ijab += (c>d) 1/16 Sum_cd gbar_cdab D_ijcd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v>v]-"), ID("[v>v]-"),
                  ID("[v,v]"), ID("[v,v]"), 1, "MO Ints <vv|vv>");
    dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "D4 <oo|vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Sigma <oo|vv>");
    dpd_contract444(&D4, &I, &S4, 0, 0, 1.0/16.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&D4);
    dpd_buf4_close(&S4);

    /*
     * S_ijab += 1/16 Sum_kl gbar_ijkl D_klab
     */
    // S_IJAB += (K>L) 1/16 Sum_KL gbar_IJKL D_KLAB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[O>O]-"),
            ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
    dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "D4 <OO|VV>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Sigma <OO|VV>");
    dpd_contract444(&I, &D4, &S4, 0, 1, 1.0/16.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&D4);
    dpd_buf4_close(&S4);

    // S_IjAb += 1/16 Sum_Kl g_IjKl D_KlAb
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D4 <Oo|Vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Sigma <Oo|Vv>");
    dpd_contract444(&I, &D4, &S4, 0, 1, 1.0/16.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&D4);
    dpd_buf4_close(&S4);

    // S_ijab += (k>l) 1/16 Sum_kl gbar_ijkl D_klab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[o>o]-"),
                  ID("[o,o]"), ID("[o,o]"), 1, "MO Ints <oo|oo>");
    dpd_buf4_init(&D4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "D4 <oo|vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Sigma <oo|vv>");
    dpd_contract444(&I, &D4, &S4, 0, 1, 1.0/16.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&D4);
    dpd_buf4_close(&S4);

    /*
     * S_ijab -= 1/16 P(ij)P(ab) Sum_kc gbar_jckb D_ikac
     */
    dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Temp (Ov|oV)");
    dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Temp (ov|ov)");
    dpd_buf4_init(&D4aa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "D4 <OO|VV>");
    dpd_buf4_sort(&D4aa, PSIF_DCFT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "D4 (OV|OV)");
    dpd_buf4_close(&D4aa);
    dpd_buf4_init(&D4aa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "D4 (OV|OV)");
    dpd_buf4_init(&D4ab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D4 <Oo|Vv>");
    dpd_buf4_sort(&D4ab, PSIF_DCFT_DPD, psqr, ID("[O,v]"), ID("[o,V]"), "D4 (Ov|oV)");
    dpd_buf4_close(&D4ab);
    dpd_buf4_init(&D4ab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "D4 (Ov|oV)");

    dpd_buf4_init(&D4bb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "D4 <oo|vv>");
    dpd_buf4_sort(&D4bb, PSIF_DCFT_DPD, prqs, ID("[o,v]"),ID("[o,v]"), "D4 (ov|ov)");
    dpd_buf4_close(&D4bb);
    dpd_buf4_init(&D4bb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "D4 (ov|ov)");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "MO Ints <oV|oV>");

    // T_IbjA = -Sum_kC D_IbkC g_jAkC
    dpd_contract444(&D4ab, &I, &Tab, 0, 0, -1.0, 0.0);

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>") ;

    // T_IbjA -= Sum_Kc g_IbKc D_KcjA
    dpd_contract444(&I, &D4ab, &Tab, 0, 1, -1.0, 1.0);

    // T_IbjA -> T_IAjb
    dpd_buf4_sort(&Tab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "Temp (OV|ov)");
    dpd_buf4_close(&Tab);
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");

    // D_IbkC -> D_ICkb
    dpd_buf4_sort(&D4ab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "D4 (OV|ov)");
    dpd_buf4_close(&D4ab);
    dpd_buf4_init(&D4ab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "D4 (OV|ov)");

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");

    // T_IAJB = Sum_kc D_IAkc g_JBkc
    dpd_contract444(&D4ab, &I, &Taa, 0, 0, 1.0, 0.0);

    // T_iajb = Sum_KC D_KCia g_KCjb
    dpd_contract444(&D4ab, &I, &Tbb, 1, 1, 1.0, 0.0);

    // T_IAjb += Sum_kc g_IAkc D_jbkc
    dpd_contract444(&I, &D4bb, &Tab, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC D_IAKC g_KCjb
    dpd_contract444(&D4aa, &I, &Tab, 0, 1, 1.0, 1.0);


    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");

    // T_IAJB -= Sum_KC D_IAKC g_JBKC
    dpd_contract444(&D4aa, &I, &Taa, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC g_IAKC D_KCjb
    dpd_contract444(&I, &D4ab, &Tab, 0, 1, -1.0, 1.0);

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");

    // T_IAJB += Sum_KC D_IAKC (JB|KC)
    dpd_contract444(&D4aa, &I, &Taa, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC (JB|KC) D_KCjb
    dpd_contract444(&I, &D4ab, &Tab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");

    // T_iajb -= Sum_kc D_iakc g_jbkc
    dpd_contract444(&D4bb, &I, &Tbb, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC D_IAkc g_jbkc
    dpd_contract444(&D4ab, &I, &Tab, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");

    // T_iajb += Sum_kc D_iakc (JB|KC)
    dpd_contract444(&D4bb, &I, &Tbb, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC D_IAkc (kc|jb)
    dpd_contract444(&D4ab, &I, &Tab, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);

    dpd_buf4_close(&D4aa);
    dpd_buf4_close(&D4ab);
    dpd_buf4_close(&D4bb);

    // T_IAJB -> T_IJAB
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    dpd_buf4_close(&Taa);
    // S_IJAB += 1/16 T_IJAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Sigma <OO|VV>");
    dpd_buf4_add(&S4, &T, 1.0/16.0);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&T);

    dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");

    // T_IJAB -> T_JIAB
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");


    // S_IJAB -= 1/16 T_JIAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Sigma <OO|VV>");
    dpd_buf4_add(&S4, &T, -1.0/16.0);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&T);

    // T_IJAB -> T_IJBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // S_IJAB -= 1/16 T_IJBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Sigma <OO|VV>");
    dpd_buf4_add(&S4, &T, -1.0/16.0);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&T);

    // T_IJAB -> T_JIBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // S_IJAB += 1/16 T_JIBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Sigma <OO|VV>");
    dpd_buf4_add(&S4, &T, 1.0/16.0);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&T);

    dpd_buf4_close(&Taa);

    // T_IAjb -> T_IjAb
    dpd_buf4_sort(&Tab, PSIF_DCFT_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    dpd_buf4_close(&Tab);
    // S_IjAb += 1/16 T_IjAb
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Sigma <Oo|Vv>");
    dpd_buf4_add(&S4, &T, 1.0/16.0);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&T);

    // T_iajb -> T_ijab
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    dpd_buf4_close(&Tbb);
    // S_ijab += 1/16 T_ijab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Sigma <oo|vv>");
    dpd_buf4_add(&S4, &T, 1.0/16.0);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&T);

    dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");

    // T_ijab -> T_jiab
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // S_ijab -= 1/16 T_jiab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Sigma <oo|vv>");
    dpd_buf4_add(&S4, &T, -1.0/16.0);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&T);

    // T_ijab -> T_ijba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // S_ijab -= 1/16 T_ijba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Sigma <oo|vv>");
    dpd_buf4_add(&S4, &T, -1.0/16.0);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&T);

    // T_ijab -> T_jiba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // S_ijab += 1/16 T_jiba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Sigma <oo|vv>");
    dpd_buf4_add(&S4, &T, 1.0/16.0);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&T);

    dpd_buf4_close(&Tbb);

}

void
DCFTSolver::compute_sigma_vector_cum_orb() {

    dpdfile2 D2, DI, DIsym;
    dpdbuf4 I, L, S4, T;

    // Lambda-Orbital block of the Hessian

    // Sigma_IJAB += 0.5 * (<CM||KI> + <CI||KM>) * Lambda_MJAB * D_KC
    // Sigma_IJAB -= 0.5 * (<CM||KJ> + <CJ||KM>) * Lambda_MIAB * D_KC

    // Sigma_IJAB -= 0.5 * (<KA||CE> + <KE||CA>) * Lambda_IJEB * D_KC;
    // Sigma_IJAB += 0.5 * (<KB||CE> + <KE||CB>) * Lambda_IJEA * D_KC;

    // Compute D * I intermediates for both terms

    // Compute D_kb * <ki||bj> intermediate (DI_ij)

    // Alpha spin

    // DI_IJ = (IJ|KC) D_KC
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "DI <O|O>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints (OO|OV)");
    dpd_contract422(&I, &D2, &DI, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // DI_IJ -= <JI|KC> D_KC
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "DI <O|O>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    dpd_contract422(&I, &D2, &DI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // DI_IJ += (IJ|ck) D_ck
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "DI <O|O>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,o]"),
                  ID("[O>=O]+"), ID("[v,o]"), 0, "MO Ints (OO|vo)");
    dpd_contract422(&I, &D2, &DI, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // Beta spin

    // DI_ij = (ij|kc) D_kc
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "DI <o|o>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints (oo|ov)");
    dpd_contract422(&I, &D2, &DI, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // DI_ij -= <ji|kc> D_kc
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "DI <o|o>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo|ov>");
    dpd_contract422(&I, &D2, &DI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // DI_ij += (ij|CK) D_CK
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "DI <o|o>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[V,O]"),
                  ID("[o,o]"), ID("[V,O]"), 0, "MO Ints (oo|VO)");
    dpd_contract422(&I, &D2, &DI, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // Compute D_kc * <ka||cb> intermediate (DI_ab)

    // Alpha spin

    // DI_AB = (AB|KC) D_KC
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "DI <V|V>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V>=V]+"), ID("[O,V]"), 0, "MO Ints (VV|OV)");
    dpd_contract422(&I, &D2, &DI, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // DI_AB -= <BA|KC> D_KC
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "DI <V|V>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V,V]"), ID("[O,V]"), 0, "MO Ints <VV|OV>");
    dpd_contract422(&I, &D2, &DI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // DI_AB += (AB|kc) D_kc
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "DI <V|V>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,v]"),
                  ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
    dpd_contract422(&I, &D2, &DI, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // Beta spin

    // DI_ab = (ab|kc) D_kc
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "DI <v|v>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,v]"),
                  ID("[v>=v]+"), ID("[o,v]"), 0, "MO Ints (vv|ov)");
    dpd_contract422(&I, &D2, &DI, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // DI_ab -= <ba|kc> D_kc
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "DI <v|v>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "D2 <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,v]"),
                  ID("[v,v]"), ID("[o,v]"), 0, "MO Ints <vv|ov>");
    dpd_contract422(&I, &D2, &DI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // DI_ab += (ab|KC) D_KC
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "DI <v|v>");
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[O,V]"),
                  ID("[v,v]"), ID("[O,V]"), 0, "MO Ints (vv|OV)");
    dpd_contract422(&I, &D2, &DI, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&D2);
    dpd_file2_close(&DI);

    // Symmetrize DI_ij and DI_ab

    // DI <O|O>
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "DI <O|O>");
    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "DI <O|O> sym");
    dpd_file2_mat_init(&DI);
    dpd_file2_mat_init(&DIsym);
    dpd_file2_mat_rd(&DI);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                DIsym.matrix[h][i][j] = DIsym.matrix[h][j][i] = DI.matrix[h][i][j] + DI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&DIsym);
    dpd_file2_close(&DI);
    dpd_file2_close(&DIsym);

    // DI <o|o>
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "DI <o|o>");
    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "DI <o|o> sym");
    dpd_file2_mat_init(&DI);
    dpd_file2_mat_init(&DIsym);
    dpd_file2_mat_rd(&DI);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                DIsym.matrix[h][i][j] = DIsym.matrix[h][j][i] = DI.matrix[h][i][j] + DI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&DIsym);
    dpd_file2_close(&DI);
    dpd_file2_close(&DIsym);

    // DI <V|V>
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "DI <V|V>");
    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "DI <V|V> sym");
    dpd_file2_mat_init(&DI);
    dpd_file2_mat_init(&DIsym);
    dpd_file2_mat_rd(&DI);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < navirpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                DIsym.matrix[h][i][j] = DIsym.matrix[h][j][i] = DI.matrix[h][i][j] + DI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&DIsym);
    dpd_file2_close(&DI);
    dpd_file2_close(&DIsym);

    // DI <v|v>
    dpd_file2_init(&DI, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "DI <v|v>");
    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "DI <v|v> sym");
    dpd_file2_mat_init(&DI);
    dpd_file2_mat_init(&DIsym);
    dpd_file2_mat_rd(&DI);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < nbvirpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                DIsym.matrix[h][i][j] = DIsym.matrix[h][j][i] = DI.matrix[h][i][j] + DI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&DIsym);
    dpd_file2_close(&DI);
    dpd_file2_close(&DIsym);

    // Sigma_IJAB += 0.5 * (<CM||KI> + <CI||KM>) * Lambda_MJAB * D_KC
    // Sigma_IJAB -= 0.5 * (<CM||KJ> + <CJ||KM>) * Lambda_MIAB * D_KC

    // Sigma_IJAB -= 0.5 * (<KA||CE> + <KE||CA>) * Lambda_IJEB * D_KC;
    // Sigma_IJAB += 0.5 * (<KB||CE> + <KE||CB>) * Lambda_IJEA * D_KC;

    //
    // Sigma <OO|VV>
    //

    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "DI <V|V> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Sigma <OO|VV>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    // The first pair index of the result of dpd_contract244 contains the file2 index.
    // The second pair index of the result of dpd_contract424 contains the file2 index.
    // This is how dpd_contract244 works: AD * IJDB -> ABIJ ?-(t)-> IJAB. If BD * IJDA -> BAIJ -(t)-> IJBA
    // This is how dpd_contract424 works: IJAD * BD -> IJAB ?-(t)-> ABIJ. If IJBD * AD -> IJBA ?-(t)-> BAIJ
    // Temp_IJAB = DIsym_AC lambda_IJCB
    dpd_contract244(&DIsym, &L, &T, 1, 2, 1, 1.0, 0.0);
    // Sigma_IJAB -= 0.5 * Temp_IJAB
    dpd_buf4_add(&S4, &T, -0.5);
    dpd_file2_close(&DIsym);
    dpd_buf4_close(&L);
    // Temp_IJAB -> Temp_IJBA
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    // Sigma_IJAB += 0.5 * Temp_IJBA
    dpd_buf4_add(&S4, &T, 0.5);
    dpd_buf4_close(&T);

    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "DI <O|O> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    // Temp_IJAB = DIsym_IK lambda_KJAB
    dpd_contract244(&DIsym, &L, &T, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&DIsym);
    dpd_buf4_close(&L);
    // Sigma_IJAB += 0.5 * Temp_IJAB
    dpd_buf4_add(&S4, &T, 0.5);
    // Temp_IJAB -> Temp_JIAB
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    // Sigma_IJAB -= 0.5 * Temp_JIAB
    dpd_buf4_add(&S4, &T, -0.5);
    dpd_buf4_close(&T);
    dpd_buf4_close(&S4);

    //
    // Sigma <Oo|Vv>
    //

    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Sigma <Oo|Vv>");
    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "DI <V|V> sym");
    // Sigma_IjAb -= 0.5 * DIsym_AC * lambda_IjCb
    dpd_contract244(&DIsym, &L, &S4, 1, 2, 1, -0.5, 1.0);
    dpd_file2_close(&DIsym);
    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "DI <v|v> sym");
    // Sigma_IjAb -= 0.5 * lambda_IjAc DIsym_bc
    dpd_contract424(&L, &DIsym, &S4, 3, 1, 0, -0.5, 1.0);
    dpd_file2_close(&DIsym);
    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "DI <O|O> sym");
    // Sigma_IjAb += 0.5 * DIsym_IK * lambda_KjAb
    dpd_contract244(&DIsym, &L, &S4, 1, 0, 0, 0.5, 1.0);
    dpd_file2_close(&DIsym);
    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "DI <o|o> sym");
    // Sigma_IjAb += 0.5 * lambda_IkAb * DIsym_jk
    dpd_contract424(&L, &DIsym, &S4, 1, 1, 1, 0.5, 1.0);
    dpd_file2_close(&DIsym);
    dpd_buf4_close(&S4);
    dpd_buf4_close(&L);

    //
    // Sigma <oo|vv>
    //

    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "DI <v|v> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Sigma <oo|vv>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    // Temp_ijab = DIsym_ac lambda_ijcb
    dpd_contract244(&DIsym, &L, &T, 1, 2, 1, 1.0, 0.0);
    // Sigma_ijab -= 0.5 * Temp_ijab
    dpd_buf4_add(&S4, &T, -0.5);
    dpd_file2_close(&DIsym);
    dpd_buf4_close(&L);
    // Temp_ijab -> Temp_ijba
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    // Sigma_ijba += 0.5 * Temp_ijba
    dpd_buf4_add(&S4, &T, 0.5);
    dpd_buf4_close(&T);

    dpd_file2_init(&DIsym, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "DI <o|o> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    // Temp_ijab = DIsym_ik lambda_kjab
    dpd_contract244(&DIsym, &L, &T, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&DIsym);
    dpd_buf4_close(&L);
    // Sigma_ijab += 0.5 * Temp_ijab
    dpd_buf4_add(&S4, &T, 0.5);
    // Temp_ijab -> Temp_jiab
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    // Sigma_ijab -= 0.5 * Temp_jiab
    dpd_buf4_add(&S4, &T, -0.5);
    dpd_buf4_close(&T);
    dpd_buf4_close(&S4);

}

int
DCFTSolver::iterate_nr_conjugate_gradients() {

    // Conjugate gradients solution of the NR equations

    double delta_new = 0.0;
    double delta_old = 0.0;
    bool converged = false;

    int cycle = 0;

    // Compute residual and form guess for the step vector D
    for (int p = 0; p < nidp_; ++p) {
        double value_r = gradient_->get(p) - sigma_->get(p) - Hd_->get(p) * X_->get(p);
        double value_d = value_r / Hd_->get(p);
        R_->set(p, value_r);
        D_->set(p, value_d);
        delta_new += value_r * value_d;
    }

    double residual_rms;

    while (!converged) {

        cycle++;

        residual_rms = 0.0;

        // Compute sigma vector
        compute_sigma_vector();

        // Compute the element of the Krylov subspace Q = Hd * D
        double dT_q = 0.0;
        for (int p = 0; p < nidp_; ++p) {
            double value_d = D_->get(p);
            Q_->set(p, sigma_->get(p) + Hd_->get(p) * value_d);
            dT_q += value_d * Q_->get(p);
        }

        // Compute scale factor for the optimal directon vector length (step length)
        double alpha = delta_new / dT_q;

        delta_old = delta_new;
        delta_new = 0.0;
        for (int p = 0; p < nidp_; ++p) {
            // Update X
            X_->set(p, X_->get(p) + alpha * D_->get(p));
            // Update the residual
            double value_r = R_->get(p) - alpha * Q_->get(p);
            R_->set(p, value_r);
            // Scale the residual to get the search direction
            double value_s = value_r / Hd_->get(p);
            S_->set(p, value_s);
            delta_new += value_r * value_s;
            residual_rms += value_r * value_r;
        }

        // Compute parameter beta for Shmidt orthogonalization
        double beta = delta_new / delta_old;

        // Compute new conjugate direction vector orthogonal to the previous search direction
        D_->scale(beta);
        D_->add(S_);

        // Compute RMS of the residual
        residual_rms = sqrt(residual_rms/nidp_);
        // Check convergence
        converged = (residual_rms < cumulant_threshold_);

        if (print_ > 3) fprintf(outfile, "%d RMS = %8.5e\n", cycle, residual_rms);
        if (cycle > maxiter_) throw PSIEXCEPTION ("Solution of the Newton-Raphson equations did not converge");

    }

    // Some attempt to use augmented Hessian approach. Think if it's needed
    double *gp = gradient_->pointer();
    double *xp = X_->pointer();
    double epsilon = C_DDOT(nidp_, gp, 1, xp, 1);
    for (int p = 0; p < nidp_; ++p) {
        X_->add(p, -epsilon * X_->get(p));
    }

    return cycle;

}

int
DCFTSolver::iterate_nr_jacobi() {

    SharedVector Xold(new Vector("Old step vector in the IDP basis", nidp_));

    bool converged_micro = false;
    int counter = 0;

    double residual_rms;

    // Jacobi solution of the NR equations
    while (!converged_micro) {

        counter++;

        residual_rms = 0.0;

        // Compute sigma vector
        compute_sigma_vector();

        double residual_rms = 0.0;
        // Update X
        for (int p = 0; p < nidp_; ++p) {
            // Compute residual: R = sigma - g + Hd * Xold
            double value_r = (-1.0) * (gradient_->get(p) - sigma_->get(p) - Hd_->get(p) * X_->get(p));
            R_->set(p, value_r);
            // Update X: Xnew = Xold - R / Hd
            if (p < orbital_idp_) {
                X_->set(p, Xold->get(p) - value_r / Hd_->get(p));
            }
            else {
                X_->set(p, Xold->get(p) - 0.25 * value_r / Hd_->get(p));
            }
            // Store the square of the residual
            residual_rms += value_r * value_r;
        }
        // Compute RMS of the residual
        residual_rms = sqrt(residual_rms/nidp_);
        // Save current X
        for (int p = 0; p < nidp_; ++p) {
            double value = X_->get(p);
            Xold->set(p, value);
            D_->set(p, value);
        }
        // Check convergence
        converged_micro = (residual_rms < cumulant_threshold_);
        if (print_ > 3) fprintf(outfile, "%d RMS = %8.5e \n", counter, residual_rms);
        if (counter > maxiter_) throw PSIEXCEPTION ("Solution of the Newton-Raphson equations did not converge");
    }

    return counter;

}

void
DCFTSolver::check_qc_convergence() {

    orbitals_convergence_ = 0.0;

    if (orbital_idp_ != 0) {
        for (int p = 0; p < orbital_idp_; ++p) orbitals_convergence_ += X_->get(p) * X_->get(p);
        orbitals_convergence_ = sqrt(orbitals_convergence_/orbital_idp_);
    }

    if(options_.get_str("QC_TYPE") == "SIMULTANEOUS") {
        cumulant_convergence_ = 0.0;

        if (cumulant_idp_ != 0) {
            for (int p = orbital_idp_; p < nidp_; ++p) cumulant_convergence_ += X_->get(p) * X_->get(p);
            cumulant_convergence_ = sqrt(cumulant_convergence_/cumulant_idp_);
        }
    }

}

void
DCFTSolver::compute_orbital_rotation_nr() {

    // Fill up the X matrix
    int orbitals_address = 0;
    int idpcount = 0;
    // Alpha spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                if (lookup_orbitals_[orbitals_address]) {
                    double value = X_->get(idpcount);
                    X_a_->set(h, i, a + naoccpi_[h], value);
                    X_a_->set(h, a + naoccpi_[h], i, (-1.0) * value);
                    idpcount++;
                }
                orbitals_address++;
            }
        }
    }

    // Beta spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = 0; a < nbvirpi_[h]; ++a){
                if (lookup_orbitals_[orbitals_address]) {
                    double value = X_->get(idpcount);
                    X_b_->set(h, i, a + nboccpi_[h], value);
                    X_b_->set(h, a + nboccpi_[h], i, (-1.0) * value);
                    idpcount++;
                }
                orbitals_address++;
            }
        }
    }

    Xtotal_a_->add(X_a_);
    Xtotal_b_->add(X_b_);

}

void
DCFTSolver::update_cumulant_nr() {

    dpdbuf4 L;

    int cumulant_address = 0;
    int idpcount = orbital_idp_;

    // Update the density cumulant
    // Alpha-Alpha spin
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&L, h);
        dpd_buf4_mat_irrep_rd(&L, h);
        #pragma omp parallel for
        for(size_t ij = 0; ij < L.params->rowtot[h]; ++ij){
            for(size_t ab = 0; ab < L.params->coltot[h]; ++ab){
                if (lookup_cumulant_[cumulant_address]) {
                    L.matrix[h][ij][ab] += 0.25 * X_->get(idpcount);
                    idpcount++;
                }
                cumulant_address++;
            }
        }
        dpd_buf4_mat_irrep_wrt(&L, h);
        dpd_buf4_mat_irrep_close(&L, h);
    }
    dpd_buf4_close(&L);

    // Alpha-Beta spin
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&L, h);
        dpd_buf4_mat_irrep_rd(&L, h);
        #pragma omp parallel for
        for(size_t ij = 0; ij < L.params->rowtot[h]; ++ij){
            for(size_t ab = 0; ab < L.params->coltot[h]; ++ab){
                if (lookup_cumulant_[cumulant_address]) {
                    L.matrix[h][ij][ab] += 0.25 * X_->get(idpcount);
                    idpcount++;
                }
                cumulant_address++;
            }
        }
        dpd_buf4_mat_irrep_wrt(&L, h);
        dpd_buf4_mat_irrep_close(&L, h);
    }
    dpd_buf4_close(&L);

    // Beta-Beta spin
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&L, h);
        dpd_buf4_mat_irrep_rd(&L, h);
        #pragma omp parallel for
        for(size_t ij = 0; ij < L.params->rowtot[h]; ++ij){
            for(size_t ab = 0; ab < L.params->coltot[h]; ++ab){
                if (lookup_cumulant_[cumulant_address]) {
                    L.matrix[h][ij][ab] += 0.25 * X_->get(idpcount);
                    idpcount++;
                }
                cumulant_address++;
            }
        }
        dpd_buf4_mat_irrep_wrt(&L, h);
        dpd_buf4_mat_irrep_close(&L, h);
    }
    dpd_buf4_close(&L);

}

void
DCFTSolver::run_davidson() {

    // Allocate the globals
    b_ = SharedMatrix(new Matrix("Expansion subspace |b>", 0, nidp_));
    vec_add_tol_ = options_.get_double("STABILITY_AUGMENT_SPACE_TOL");
    r_convergence_ = options_.get_double("STABILITY_CONVERGENCE");
    n_add_ = options_.get_int("STABILITY_ADD_VECTORS");
    nguess_ = options_.get_int("STABILITY_N_GUESS_VECTORS");
    nevals_ = options_.get_int("STABILITY_N_EIGENVALUES");
    max_space_ = options_.get_int("STABILITY_MAX_SPACE_SIZE");
    b_dim_ = 0;

    // Create a set of guess orthonormal expansion vectors b_ (identity matrix)
    davidson_guess();
    bool converged = false;
    int count = 0;

    // Create a matrix with the Hessian diagonal for convenience
    fprintf(outfile, "\tStability analysis of DCFT solution \n");
    SharedVector Evals;
    SharedMatrix Evecs;
    while(!converged){
        if (count > maxiter_) throw PSIEXCEPTION("Davidson diagonalization did not converge!");
        if (print_ > 1) fprintf(outfile, "\tIteration %d\n", ++count);

        // Form the off-diagonal contribution to the sigma vector
        SharedMatrix sigma_vector(new Matrix("Sigma vector for the Davidson algorithm", b_dim_, nidp_));
        for (int k = 0; k < b_dim_; ++k) {
            double *D_p = D_->pointer();
            double *b_p = b_->pointer()[k];
            // Copy the kth b vector to D_
            ::memcpy(D_p, b_p, nidp_ * sizeof(double));
            // Compute the off-diagonal part of the sigma vector for kth b vector sigma = Ho * b
            compute_sigma_vector();
            // Store sigma vector
            double *sigma_p = sigma_->pointer();
            double *sigma_vector_p = sigma_vector->pointer()[k];
            ::memcpy(sigma_vector_p, sigma_p, nidp_ * sizeof(double));
        }
        // Add the diagonal contribution
        SharedMatrix temp(new Matrix("Temp matrix", b_dim_, nidp_));
        temp->copy(b_);
        double **temp_p = temp->pointer();
        for (int i = 0; i < nidp_; ++i) C_DSCAL(b_dim_, Hd_->get(i), &temp_p[0][i], nidp_);

        sigma_vector->add(temp);

        // Form G matrix and diagonalize it
        SharedMatrix G(new Matrix("Subspace representation of the Hessian", b_dim_, b_dim_));
        Evecs = SharedMatrix(new Matrix("Eigenvectors of the Hessian subspace representation", b_dim_, b_dim_));
        Evals = SharedVector(new Vector("Eigenvalues of the Hessian subspace representation", b_dim_));
        G->gemm(false, true, 1.0, b_, sigma_vector, 0.0);
        G->diagonalize(Evecs, Evals, ascending);

        // Define the eigenvectors to be positive to make sure the phase doesn't change
        double *ones = new double[b_dim_];
        for (int i = 0; i < b_dim_; ++i) ones[i] = 1.0;
        double **Evecs_p = Evecs->pointer();
        for (int k = 0; k < b_dim_; ++k) {
            double dot = C_DDOT(b_dim_, ones, 1, &Evecs_p[0][k], b_dim_);
            if (dot < 0.0) Evecs->scale_column(0, k, -1.0);
        }

        // Rotate sigma and b vectors to the new subspace
        // sigma'_kp = sum(i) alpha_ki sigma_ip
        // b_kp' = sum(i) alpha_ki b_ip
        temp->gemm(true, false, 1.0, Evecs, sigma_vector, 0.0);
        sigma_vector->copy(temp);
        temp->gemm(true, false, 1.0, Evecs, b_, 0.0);
        b_->copy(temp);

        // Check whether the norm is preserved
        SharedMatrix check(new Matrix("Orthonormality Check, after rotation", b_dim_, b_dim_));
        check->gemm(0, 1, 1.0, b_, b_, 0.0);
        for (int i = 0; i < b_dim_; ++i) {
            if (fabs(check->get(i,i) - 1.0) > 1e-5) throw PSIEXCEPTION("Norm is not preserved! Make STABILITY_AUGMENT_SPACE_TOL larger");
        }

        // Compute the residual: r_kp = sigma'_kp - ro_k b_kp'
        SharedMatrix r(new Matrix("Residual", b_dim_, nidp_));
        r->copy(sigma_vector);
        double **r_p = r->pointer();
        double **b_p = b_->pointer();
        for (int k = 0; k < b_dim_; ++k) C_DAXPY(nidp_, -Evals->get(k), b_p[k], 1, r_p[k], 1);

        // Check the convergence for each vector
        double max_rms = 0.0;
        int n_good = 0;
        int n_bad  = 0;
        int nvecs = nevals_ < b_dim_ ? nevals_ : b_dim_;
        if (print_ > 1) fprintf(outfile, "\tEigenvalues:\n");
        for (int k = 0; k < nvecs; ++k) {
            double new_val = Evals->get(k);
            double ms = C_DDOT(b_dim_, r_p[k], 1, r_p[k], 1);
            double rms = sqrt(ms / (double) b_dim_);
            bool not_converged = rms > r_convergence_;
            if (not_converged) n_bad += 1;
            else n_good += 1;
            max_rms = rms > max_rms ? rms : max_rms;
            if (print_ > 1) fprintf(outfile,  "\t\t%s%10.6f   (residual = %8.3e)\n", not_converged ? "*" : " ", new_val, rms);
        }

        if (print_ > 1) fprintf(outfile, "\tThere are %d vectors in the subspace\n", b_dim_);
        if (print_ > 1) fprintf(outfile, "\tMax RMS residual %8.3e, %d converged, %d not converged\n", max_rms, n_good, n_bad);

        if (n_bad == 0) {
            converged = true;
            break;
        }

        // Obtain new vectors from the residual: delta = -r_kp / (Hd_p - Evals_k)
        SharedMatrix delta(new Matrix ("Correction vector", b_dim_, nidp_));
        for (int k = 0; k < b_dim_; ++k)
            for (int p = 0; p < nidp_; ++p)
                delta->set(k, p, -r->get(k, p));
//                delta->set(k, p, -r->get(k, p) / (Hd_->get(p) - Evals->get(k)));

        // Orthogonalize the new vectors against the subspace, and add if nececssary
        int added_vectors = 0;
        double **delta_p = delta->pointer();
        int subspace_size = b_dim_;
        for(int k = subspace_size - 1; k >= 0; --k){
            if(augment_b(delta_p[k], vec_add_tol_)){
                added_vectors++;
                if(added_vectors == n_add_)
                    break;
            }
        }
        if(!added_vectors && !converged) throw PSIEXCEPTION("No new vectors were added");

        if (print_ > 1) fprintf(outfile, "\tAdded %d new vector(s) to the subspace\n\n\n", added_vectors);

        if (b_dim_ > max_space_) throw PSIEXCEPTION("The subspace size is exceeded, but the convergence is not reached in stability analysis");

    }

    // Perform the stability analysis by analyzing the eigenvector of the Hessian for several largest contributions
    fprintf(outfile, "\tLowest %d eigenvalues of the electronic Hessian: \n", nevals_);
    int values_to_print = 5;
    int n_neg = 0;
    int nvecs = nevals_ < b_dim_ ? nevals_ : b_dim_;
    for (int k = 0; k < nvecs; k++) {
        double value = Evals->get(k);
        fprintf(outfile, "\t %10.6f \n", value);
        if (value < 0.0) {
            double *max_values = new double[values_to_print + 1];
            int *max_values_idp = new int[values_to_print + 1];
            int stored = 0;

            double norm = 0.0;

            for (int i = 0; i < nidp_; ++i) {
                double evec_value = b_->get(k, i);
                norm += evec_value * evec_value;
                // Store the value at the end of the list
                max_values[stored] = evec_value;
                max_values_idp[stored] = i;
                // Sort the values in the ascending order
                for (int p = 0; p < stored; ++p) {
                    if (fabs(max_values[stored - p]) > fabs(max_values[stored - p - 1])) {
                        double temp_val = max_values[stored - p - 1];
                        max_values[stored - p - 1] = max_values[stored - p];
                        max_values[stored - p] = temp_val;
                        int temp_val_idp = max_values_idp[stored - p - 1];
                        max_values_idp[stored - p - 1] = max_values_idp[stored - p];
                        max_values_idp[stored - p] = temp_val_idp;
                    }
                    else break;
                }
                if (stored < values_to_print) stored++;

            }
            fprintf(outfile, "\t %d largest contributions to the eigenvector: \n", values_to_print);
            for (int i = 0; i < values_to_print; ++i) {
                fprintf(outfile, "\t %10.3e %s \n", max_values[i], (max_values_idp[i] < orbital_idp_) ? ("orbital space") : ("cumulant space"));
            }
            n_neg++;
        }
    }
    if (n_neg) fprintf(outfile, "\tSolution is unstable (%d negative eigenvalues obtained) \n", n_neg);

}

void
DCFTSolver::davidson_guess() {

    int count = 0;
    int dimension = nguess_ < nidp_ ? nguess_ : nidp_;
    while (count < dimension) {
        Vector temp("Temp", nidp_);
        temp.set(count, 1.0);
        // To avoid singularity in the solution the second element is set to be non-zero
        temp.set(count + 1, 0.1);
        double *ptemp = temp.pointer();
        // Orthonormalize the guess vectors and form the guess b
        if (augment_b(ptemp, vec_add_tol_)) count++;
    }

}

bool
DCFTSolver::augment_b(double *vec, double tol) {

    // Normalize the vec array first
    double vec_norm = sqrt(C_DDOT(nidp_, vec, 1, vec, 1));
    double inv_norm = 1.0/vec_norm;
    C_DSCAL(nidp_, inv_norm, vec, 1);

    // Allocate |b'> vector and copy vec to it
    SharedMatrix bprime(new Matrix("B'", 1, nidp_));
    double **bpp = bprime->pointer();
    ::memcpy(bpp[0], vec, nidp_ * sizeof(double));

    // Compute the overlap of the new |b'> space and the existing space |b>: <b|b'>
    SharedMatrix bxb(new Matrix("<b'|b>", b_dim_, 1));
    if (b_dim_) bxb->gemm(0, 1, 1.0, b_, bprime, 0.0);

    // Project out the existing subspace |b> from the current |b'> subspace: |b'> = |b'> - |b><b|b'>
    for (int k = 0; k < b_dim_; ++k) C_DAXPY(nidp_, -bxb->get(k, 0), b_->pointer()[k], 1, bpp[0], 1);

    // Compute the norm of the leftover part of the |b'> space
    vec_norm = sqrt(C_DDOT(nidp_, bpp[0], 1, bpp[0], 1));

    // If the rms exceeds a given threshold then augment the |b> space
    if(vec_norm > tol){
        double inv_norm = 1.0 / vec_norm;
        C_DSCAL(nidp_, inv_norm, bpp[0], 1);

        std::vector<SharedMatrix> mats;
        mats.push_back(b_);
        mats.push_back(bprime);
        b_ = Matrix::vertcat(mats);
        b_->set_name("B");
        b_dim_++;
        return true;
    }

    return false;

}

}} //End namespaces

