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

    bool orbitalsDone     = false;
    bool cumulantDone     = false;
    bool energyConverged  = false;
    bool densityConverged = false;

    int cycle = 0;

    if (options_.get_str("DERTYPE") == "FIRST") throw PSIEXCEPTION ("Analytic gradients are not available for QC DCFT yet. Please choose a different algorithm");

    // Allocate the memory
    qc_dcft_init();

    while((!orbitalsDone || !cumulantDone || !energyConverged || !densityConverged) && cycle++ < maxiter_){

        // Compute the generalized Fock matrix and orbital gradient ([F,Kappa]) in the MO basis
        compute_orbital_gradient();
        // Build G and F intermediates needed for the density cumulant residual equations and DCFT energy computation
        build_intermediates();
        // Compute the residuals for density cumulant equations
        lambda_convergence_ = compute_lambda_residual();
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
        energyConverged = fabs(old_total_energy_ - new_total_energy_) < lambda_threshold_;
        // Print the iterative trace
        fprintf(outfile, "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f       *\n",
                cycle, scf_convergence_, lambda_convergence_, new_total_energy_ - old_total_energy_,
                new_total_energy_);
        fflush(outfile);
        // Determine the independent pairs (IDPs) and create array for the orbital and cumulant gradient in the basis of IDPs
        form_idps();
        if (nidp_ == 0) {
            orbitalsDone = true; cumulantDone = true; energyConverged = true; densityConverged = true;
            break;
        }
        // Compute sigma vector in the basis of IDPs
        compute_sigma_vector();
        // Solve the NR equations using conjugate gradients
        iterate_conjugate_gradients();
        // Check the convergence by computing the change in the orbitals and the cumulant
        check_qc_convergence();
        orbitalsDone = scf_convergence_ < scf_threshold_;
        cumulantDone = lambda_convergence_ < lambda_threshold_;
        // Update cumulant and orbitals with the solution of the Newton-Raphson equations
        update_cumulant_and_orbitals();
        if (orbital_idp_ != 0) {
            // Update the density
            densityConverged = update_scf_density() < scf_threshold_;
            // Write orbitals to the checkpoint file
            write_orbitals_to_checkpoint();
            // Transform two-electron integrals to the MO basis using new orbitals, build denominators
            // TODO: Transform_integrals shouldn't call build denominators for the QC alogorithm
            transform_integrals();
        }
    }

    if(!orbitalsDone || !cumulantDone || !densityConverged)
        throw ConvergenceError<int>("DCFT", maxiter_, lambda_threshold_,
                               lambda_convergence_, __FILE__, __LINE__);

    // Make sure that the orbital phase is retained and the Fock matrix is diagonal for the gradients

}

void
DCFTSolver::qc_dcft_init(){

    orbital_gradient_a_ = SharedMatrix(new Matrix("MO basis Orbital Gradient (Alpha)", nirrep_, nmopi_, nmopi_));
    orbital_gradient_b_ = SharedMatrix(new Matrix("MO basis Orbital Gradient (Beta)", nirrep_, nmopi_, nmopi_));

    // The number of IDPs is set to zero in the beginning
    nidp_ = 0;

    dim_ = nalpha_ * navir_ + nbeta_ * nbvir_;
    dim_ += (nalpha_ * (nalpha_ - 1) / 2) * (navir_ * (navir_ - 1) / 2);
    dim_ += (nalpha_ * nbeta_) * (navir_ * nbvir_);
    dim_ += (nbeta_ * (nbeta_ - 1) / 2) * (nbvir_ * (nbvir_ - 1) / 2);

    lookup_ = new int[dim_];
    ::memset(lookup_, '\0', sizeof(int)*dim_);

}

void
DCFTSolver::compute_orbital_gradient(){

    // Initialize the idempotent contribution to the OPDM (Kappa)
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
    // Build guess Tau from the density cumulant in the MO basis and transform it to the SO basis
    build_tau();
    transform_tau();
    // Copy core hamiltonian into the Fock matrix array: F = H
    Fa_->copy(so_h_);
    Fb_->copy(so_h_);
    // Build the new Fock matrix from the SO integrals: F += Gbar * Kappa
    process_so_ints();
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

    // Form F * kappa (Alpha spin)
    orbital_gradient_a_->gemm(false, false, 2.0, moFa_, full_kappa_a, 0.0);
    // Form -kappa * F (Alpha spin) to obtain [F,kappa]
    orbital_gradient_a_->gemm(false, false, -2.0, full_kappa_a, moFa_, 1.0);

    // Form F * kappa (Beta spin)
    orbital_gradient_b_->gemm(false, false, 2.0, moFb_, full_kappa_b, 0.0);
    // Form -kappa * F (Alpha spin) to obtain [F,kappa]
    orbital_gradient_b_->gemm(false, false, -2.0, full_kappa_b, moFb_, 1.0);

}

void DCFTSolver::form_idps(){

    // Ignore orbital gradient elements that are less than this value
    double cutoff = ((1.0e-10 < (lambda_threshold_  * 0.01)) ? 1.0e-10 : (lambda_threshold_  * 0.01));

    // Zero out the counters
    int old_nidp = nidp_;
    nidp_ = 0;
    orbital_idp_a_ = 0;
    orbital_idp_b_ = 0;
    orbital_idp_ = 0;
    lambda_idp_aa_ = 0;
    lambda_idp_ab_ = 0;
    lambda_idp_bb_ = 0;
    lambda_idp_ = 0;

    // Zero the lookup array
    ::memset(lookup_, '\0', sizeof(int)*dim_);

    // Temporary vectors containing gradient value and diagonal part of the Hessian for each IDP
    double *grad = new double[dim_];
    double *Hd = new double[dim_];

    // Count the number of IDPs for orbital rotations (Alpha spin)
    // The minus sign in the gradient takes into account the sign of the g vector in the N-R equations: dX H = -g
    int address = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                if (fabs(orbital_gradient_a_->get(h, i, a + naoccpi_[h])) > cutoff) {
                    lookup_[address] = 1;
                    grad[orbital_idp_a_] = (-1.0) * orbital_gradient_a_->get(h, i, a + naoccpi_[h]);
                    Hd[orbital_idp_a_] = 2.0 * (moFa_->get(h, a + naoccpi_[h], a + naoccpi_[h]) - moFa_->get(h, i, i));
                    orbital_idp_a_++;
                }
                address++;
            }
        }
    }

    // Count the number of IDPs for orbital rotations (Beta spin)
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = 0; a < nbvirpi_[h]; ++a){
                if (fabs(orbital_gradient_b_->get(h, i, a + nboccpi_[h])) > cutoff) {
                    lookup_[address] = 1;
                    int index = orbital_idp_a_ + orbital_idp_b_;
                    grad[index] = (-1.0) * orbital_gradient_b_->get(h, i, a + nboccpi_[h]);
                    Hd[index] = 2.0 * (moFb_->get(h, a + nboccpi_[h], a + nboccpi_[h]) - moFb_->get(h, i, i));
                    orbital_idp_b_++;
                }
                address++;
            }
        }
    }

    orbital_idp_ = orbital_idp_a_ + orbital_idp_b_;

    // Count the number of IDPs for cumulant updates

    dpdbuf4 R, I;

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
                    lookup_[address] = 1;
                    int index = orbital_idp_ + lambda_idp_aa_;
                    grad[index] = -0.25 * R.matrix[h][ij][ab];
                    double value = moFa_->get(asym, a + naoccpi_[asym], a + naoccpi_[asym])
                            + moFa_->get(bsym, b + naoccpi_[bsym], b + naoccpi_[bsym])
                            - moFa_->get(isym, i, i) - moFa_->get(jsym, j, j);
                    Hd[index] = (1.0/16.0) * value;
                    lambda_idp_aa_++ ;
                }
                address++;
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
                    lookup_[address] = 1;
                    int index = orbital_idp_ + lambda_idp_aa_ + lambda_idp_ab_;
                    grad[index] = -0.25 * R.matrix[h][ij][ab];
                    double value = moFa_->get(asym, a + naoccpi_[asym], a + naoccpi_[asym])
                            + moFb_->get(bsym, b + nboccpi_[bsym], b + nboccpi_[bsym])
                            - moFa_->get(isym, i, i) - moFb_->get(jsym, j, j);
                    Hd[index] = (1.0/16.0) * value;
                    lambda_idp_ab_++ ;
                }
                address++;
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
                    lookup_[address] = 1;
                    int index = orbital_idp_ + lambda_idp_aa_ + lambda_idp_ab_ + lambda_idp_bb_;
                    grad[index] = -0.25 * R.matrix[h][ij][ab];
                    double value = moFb_->get(asym, a + nboccpi_[asym], a + nboccpi_[asym])
                            + moFb_->get(bsym, b + nboccpi_[bsym], b + nboccpi_[bsym])
                            - moFb_->get(isym, i, i) - moFb_->get(jsym, j, j);
                    Hd[index] = (1.0/16.0) * value;
                    lambda_idp_bb_++ ;
                }
                address++;
            }
        }
        dpd_buf4_mat_irrep_close(&R, h);
    }
    dpd_buf4_close(&R);

    lambda_idp_ = lambda_idp_aa_ + lambda_idp_ab_ + lambda_idp_bb_;
    nidp_ = orbital_idp_ + lambda_idp_;

    // Reallocate the memory for the arrays of the dimensions changed
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
    sigma_->zero();

    // Set the gradient Hessian diagonal vectors
    for (int p = 0; p < nidp_; ++p) {
        gradient_->set(p, grad[p]);
        Hd_->set(p, Hd[p]);
        // Compute the guess for the step vector
        X_->set(p, grad[p]/Hd[p]);
        // D is used to store X for the initial formation of the sigma vector
        D_->set(p, grad[p]/Hd[p]);
    }

}

void
DCFTSolver::compute_sigma_vector() {

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdfile2 D2, S2, LD_OO, LD_oo, LD_VV, LD_vv, DI, DIsym;
    dpdbuf4 I, L, D4, S4, T, Taa, Tab, Tbb, D4aa, D4ab, D4bb;

    // Copy the IDP arrays from memory to DPD file2 and buf4

    int address = 0;
    int idpcount = 0;
    // Conjugate direction for orbital rotations (Alpha spin)
    dpd_file2_init(&D2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "D2 <O|V>");
    dpd_file2_mat_init(&D2);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                D2.matrix[h][i][a] = 0.0;
                if (lookup_[address]) {
                    D2.matrix[h][i][a] = D_->get(idpcount);
                    idpcount++;
                }
                address++;
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
                if (lookup_[address]) {
                    D2.matrix[h][i][a] = D_->get(idpcount);
                    idpcount++;
                }
                address++;
            }
        }
    }
    dpd_file2_mat_wrt(&D2);
    dpd_file2_close(&D2);

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
                if (lookup_[address]) {
                    D4.matrix[h][ij][ab] = D_->get(idpcount);
                    idpcount++;
                }
                address++;
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
                if (lookup_[address]) {
                    D4.matrix[h][ij][ab] = D_->get(idpcount);
                    idpcount++;
                }
                address++;
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
                if (lookup_[address]) {
                    D4.matrix[h][ij][ab] = D_->get(idpcount);
                    idpcount++;
                }
                address++;
            }
        }
        dpd_buf4_mat_irrep_wrt(&D4, h);
        dpd_buf4_mat_irrep_close(&D4, h);
    }
    dpd_buf4_close(&D4);

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

    if (options_.get_bool("QC_COUPLING")) {
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

    if (options_.get_bool("QC_COUPLING")) {
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

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    //
    // Now copy the sigma vector back to memory
    //

    address = 0;
    idpcount = 0;
    // Alpha spin
    dpd_file2_init(&S2, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Sigma <O|V>");
    dpd_file2_mat_init(&S2);
    dpd_file2_mat_rd(&S2);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                if (lookup_[address]) {
                    sigma_->set(idpcount, S2.matrix[h][i][a]);
                    idpcount++;
                }
                address++;
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
                if (lookup_[address]) {
                    sigma_->set(idpcount, S2.matrix[h][i][a]);
                    idpcount++;
                }
                address++;
            }
        }
    }
    dpd_file2_close(&S2);

    // Alpha-Alpha spin
    dpd_buf4_init(&S4, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
              ID("[O>O]-"), ID("[V>V]-"), 0, "Sigma <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&S4, h);
        dpd_buf4_mat_irrep_rd(&S4, h);
        #pragma omp parallel for
        for(size_t ij = 0; ij < S4.params->rowtot[h]; ++ij){
            for(size_t ab = 0; ab < S4.params->coltot[h]; ++ab){
                if (lookup_[address]) {
                    sigma_->set(idpcount, S4.matrix[h][ij][ab]);
                    idpcount++;
                }
                address++;
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
                if (lookup_[address]) {
                    sigma_->set(idpcount, S4.matrix[h][ij][ab]);
                    idpcount++;
                }
                address++;
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
                if (lookup_[address]) {
                    sigma_->set(idpcount, S4.matrix[h][ij][ab]);
                    idpcount++;
                }
                address++;
            }
        }
        dpd_buf4_mat_irrep_close(&S4, h);
    }
    dpd_buf4_close(&S4);

}

void
DCFTSolver::iterate_conjugate_gradients() {

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
        converged = (residual_rms < lambda_threshold_);

        cycle++;
        fprintf(outfile, "%d RMS = %8.5e\n", cycle, residual_rms);
        if (cycle > maxiter_) throw PSIEXCEPTION ("Solution of the Newton-Raphson equations did not converge");

    }

    // Some attempt to use augmented Hessian approach. Think if it's needed
    double *gp = gradient_->pointer();
    double *xp = X_->pointer();
    double epsilon = C_DDOT(nidp_, gp, 1, xp, 1);
    for (int p = 0; p < nidp_; ++p) {
        X_->add(p, -epsilon * X_->get(p));
    }


}

void
DCFTSolver::check_qc_convergence() {

    scf_convergence_ = 0.0;
    lambda_convergence_ = 0.0;

    if (orbital_idp_ != 0) {
        for (int p = 0; p < orbital_idp_; ++p) scf_convergence_ += X_->get(p) * X_->get(p);
        scf_convergence_ = sqrt(scf_convergence_/orbital_idp_);
    }

    if (lambda_idp_ != 0) {
        for (int p = orbital_idp_; p < nidp_; ++p) lambda_convergence_ += X_->get(p) * X_->get(p);
        lambda_convergence_ = sqrt(lambda_convergence_/lambda_idp_);
    }

}

void
DCFTSolver::update_cumulant_and_orbitals() {

    dpdbuf4 L;

    // Initialize the orbital rotation matrix
    SharedMatrix U_a(new Matrix("Orbital rotation matrix (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedMatrix U_b(new Matrix("Orbital rotation matrix (Beta)", nirrep_, nmopi_, nmopi_));
    SharedMatrix X_a(new Matrix("Generator of the orbital rotations (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedMatrix X_b(new Matrix("Generator of the orbital rotations (Beta)", nirrep_, nmopi_, nmopi_));

    // Fill up the X matrix
    int address = 0;
    int idpcount = 0;
    // Alpha spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                if (lookup_[address]) {
                    double value = X_->get(idpcount);
                    X_a->set(h, i, a + naoccpi_[h], value);
                    X_a->set(h, a + naoccpi_[h], i, (-1.0) * value);
                    idpcount++;
                }
                address++;
            }
        }
    }

    // Beta spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = 0; a < nbvirpi_[h]; ++a){
                if (lookup_[address]) {
                    double value = X_->get(idpcount);
                    X_b->set(h, i, a + nboccpi_[h], value);
                    X_b->set(h, a + nboccpi_[h], i, (-1.0) * value);
                    idpcount++;
                }
                address++;
            }
        }
    }

    if (orbital_idp_ != 0) {
        // U = I
        U_a->identity();
        U_b->identity();

        // U += X
        U_a->add(X_a);
        U_b->add(X_b);

        // U += 0.5 * X * X
        U_a->gemm(false, false, 0.5, X_a, X_a, 1.0);
        U_b->gemm(false, false, 0.5, X_b, X_b, 1.0);

        // Orthogonalize the U vectors
        int rowA = U_a->nrow();
        int colA = U_a->ncol();

        double **U_a_block = block_matrix(rowA, colA);
        memset(U_a_block[0], 0, sizeof(double)*rowA*colA);
        U_a_block = U_a->to_block_matrix();
        schmidt(U_a_block, rowA, colA, outfile);
        U_a->set(U_a_block);
        free_block(U_a_block);

        int rowB = U_a->nrow();
        int colB = U_b->ncol();

        double **U_b_block = block_matrix(rowB, colB);
        memset(U_b_block[0], 0, sizeof(double)*rowB*colB);
        U_b_block = U_b->to_block_matrix();
        schmidt(U_b_block, rowB, colB, outfile);
        U_b->set(U_b_block);
        free_block(U_b_block);

        // Rotate the orbitals
        old_ca_->copy(Ca_);
        old_cb_->copy(Cb_);

        Ca_->gemm(false, false, 1.0, old_ca_, U_a, 0.0);
        Cb_->gemm(false, false, 1.0, old_cb_, U_b, 0.0);
    }

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
                if (lookup_[address]) {
                    L.matrix[h][ij][ab] += 0.25 * X_->get(idpcount);
                    idpcount++;
                }
                address++;
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
                if (lookup_[address]) {
                    L.matrix[h][ij][ab] += 0.25 * X_->get(idpcount);
                    idpcount++;
                }
                address++;
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
                if (lookup_[address]) {
                    L.matrix[h][ij][ab] += 0.25 * X_->get(idpcount);
                    idpcount++;
                }
                address++;
            }
        }
        dpd_buf4_mat_irrep_wrt(&L, h);
        dpd_buf4_mat_irrep_close(&L, h);
    }
    dpd_buf4_close(&L);

}

}} //End namespaces

