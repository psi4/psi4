/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#ifndef _PSI4_SRC_BIN_DCT_DCT_H_
#define _PSI4_SRC_BIN_DCT_DCT_H_

#include "psi4/dct/dct_df_tensor.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/dimension.h"

#define DCT_TIMER

// Handy mints timer macros, requires libqt to be included
#ifdef DCT_TIMER
#include "psi4/libqt/qt.h"

#define dct_timer_on(a) timer_on((a));
#define dct_timer_off(a) timer_off((a));
#else
#define dct_timer_on(a)
#define dct_timer_off(a)
#endif

#define ID(x) _ints->DPD_ID(x)

#ifndef INDEX
#define INDEX(i, j) (((i) > (j)) ? (((i) * ((i) + 1) / 2) + (j)) : (((j) * ((j) + 1) / 2) + (i)))
#endif

#define PRINT_ENERGY_COMPONENTS 0

namespace psi {

class Options;
class IntegralTransform;

namespace dct {

class DCTSolver : public Wavefunction {
   public:
    DCTSolver(SharedWavefunction ref_wfn, Options& options);
    ~DCTSolver() override;

    double compute_energy() override;

   protected:
    std::unique_ptr<IntegralTransform> _ints;

    void compute_relaxed_opdm();
    void initialize_amplitudes();
    void initialize_orbitals_from_reference_U();
    void initialize_orbitals_from_reference_R();
    void initialize_integraltransform();
    void finalize();
    void transform_integrals();
    void transform_core_integrals();
    void sort_OOOO_integrals();
    void sort_OOVV_integrals();
    void sort_VVVV_integrals();
    void sort_OVOV_integrals();
    void sort_OVVV_integrals();
    void sort_OOOV_integrals();
    void init();
    void compute_dct_energy();
    void compute_cepa0_energy();
    void update_cumulant_jacobi();
    void compute_scf_energy();
    void compute_SO_tau_U();
    void compute_SO_tau_R();
    void build_d_fourth_order_U();
    void build_d_U();
    void build_d_R();
    void build_tau_U();
    void build_tau_R();
    void transform_tau_U();
    void transform_tau_R();
    void build_gtau();
    void print_opdm();
    void build_cumulant_intermediates();
    void process_so_ints();
    void build_AO_tensors();
    void build_denominators();
    void update_fock();
    void dump_density();
    void dpd_buf4_add(dpdbuf4* A, dpdbuf4* B, double alpha);
    void half_transform(dpdbuf4* A, dpdbuf4* B, SharedMatrix& C1, SharedMatrix& C2, int* mospi_left, int* mospi_right,
                        int** so_row, int** mo_row, bool backwards, double alpha, double beta);
    void file2_transform(dpdfile2* A, dpdfile2* B, SharedMatrix C, bool backwards);
    void AO_contribute(dpdbuf4* tau1_AO, dpdbuf4* tau2_AO, int p, int q, int r, int s, double value,
                       dpdfile2* = nullptr, dpdfile2* = nullptr, dpdfile2* = nullptr);
    // void AO_contribute(dpdfile2 *tau1_AO, dpdfile2 *tau2_AO, int p, int q,
    //        int r, int s, double value);
    bool correct_mo_phases(bool dieOnError = true);
    bool correct_mo_phase_spincase(Matrix& temp, Matrix& overlap, const Matrix& old_C, Matrix& C, bool dieOnError = true) const;
    double compute_cumulant_residual();
    double compute_scf_error_vector();
    double update_scf_density(bool damp = false);
    void run_twostep_dct();
    int run_twostep_dct_cumulant_updates();
    void run_twostep_dct_orbital_updates();
    void run_simult_dct();
    void run_simult_dct_oo();
    // DCT analytic gradient subroutines
    SharedMatrix compute_gradient() override;
    void dc06_response_init();
    void response_guess();
    void dc06_response();
    void dc06_compute_relaxed_density_1PDM();
    void oo_gradient_init();
    void compute_lagrangian_OV();
    void compute_lagrangian_VO();
    void iterate_orbital_response();
    void orbital_response_guess();
    void compute_orbital_response_intermediates();
    double update_orbital_response();
    double compute_response_coupling();
    void iterate_cumulant_response();
    void cumulant_response_guess();
    void build_perturbed_tau();
    void compute_cumulant_response_intermediates();
    double compute_cumulant_response_residual();
    void update_cumulant_response();
    void compute_lagrangian_OO(bool separate_gbargamma);
    void compute_lagrangian_VV(bool separate_gbargamma);
    void compute_ewdm_dc();
    void compute_ewdm_odc();
    void compute_relaxed_density_OOOO();
    void compute_relaxed_density_OOVV();
    void compute_relaxed_density_OVOV();
    void compute_relaxed_density_VVVV();
    void compute_TPDM_trace(bool cumulant_only);
    // Quadratically-convergent DCT
    void run_qc_dct();
    void compute_orbital_gradient();
    void form_idps();
    void compute_sigma_vector();
    void compute_sigma_vector_orb_orb();
    void compute_sigma_vector_orb_cum();
    void compute_sigma_vector_cum_cum();
    void compute_sigma_vector_cum_orb();
    int iterate_nr_conjugate_gradients();
    int iterate_nr_jacobi();
    void check_qc_convergence();
    void compute_orbital_rotation_nr();
    void update_cumulant_nr();
    void run_davidson();
    void davidson_guess();
    // Exact Tau
    void refine_tau();
    void form_density_weighted_fock();
    // Cumulant residual intermediates
    void compute_G_intermediate();
    void compute_F_intermediate();
    void compute_V_intermediate();
    void compute_W_intermediate();
    void compute_H_intermediate();
    void compute_I_intermediate();
    void compute_J_intermediate();
    void compute_K_intermediate();
    void compute_L_intermediate();
    void compute_O_intermediate();
    void compute_M_intermediate();
    void compute_N_intermediate();

    // Orbital-optimized DCT
    /// Converge a DC-12 computation to obtain guess orbitals and double amplitudes
    void run_simult_dc_guess();
    double compute_orbital_residual();
    // Compute 2RDMs. These arise as intermediates in computing the
    // orbital residual. For DF variants, only the cumulant of the
    // 2RDM is needed (as the RIFIT part).
    // Setting cumulant_only saves as "Cumulant", else "Gamma", in
    // addition to changing what is actually computed.
    void compute_unrelaxed_density_OOOO(bool cumulant_only);
    void compute_unrelaxed_separable_density_OOOO();
    void compute_unrelaxed_density_OOVV(bool cumulant_only);
    void compute_unrelaxed_density_OVOV(bool cumulant_only);
    void compute_unrelaxed_separable_density_OVOV();
    void compute_unrelaxed_density_VVVV(bool cumulant_only);
    void compute_unrelaxed_separable_density_VVVV();
    void compute_orbital_gradient_OV(bool separate_gbargamma);
    void compute_orbital_gradient_VO(bool separate_gbargamma);
    void compute_orbital_rotation_jacobi();
    /// target = old * exp(X)
    void rotate_matrix(const Matrix& X, const Matrix& old, Matrix& target);
    void rotate_orbitals();
    Matrix construct_oo_density(const Matrix& occtau, const Matrix& virtau, const Matrix& kappa, const Matrix& C);
    // Three-particle cumulant contributions
    double compute_three_particle_energy();
    void dct_semicanonicalize();
    void transform_integrals_triples();
    void dump_semicanonical();
    void semicanonicalize_gbar_ovvv();
    void semicanonicalize_gbar_ooov();
    void semicanonicalize_dc();
    double compute_triples_aaa();
    double compute_triples_aab();
    double compute_triples_abb();
    double compute_triples_bbb();

    // RHF-reference DCT
    double compute_energy_RHF();
    double update_scf_density_RHF(bool damp = false);
    double compute_scf_error_vector_RHF();
    void build_denominators_RHF();
    void initialize_amplitudes_RHF();
    void transform_integrals_RHF();
    void transform_core_integrals_RHF();
    void sort_OOOO_integrals_RHF();
    void sort_OOVV_integrals_RHF();
    void sort_VVVV_integrals_RHF();
    void sort_OVOV_integrals_RHF();
    void sort_OVVV_integrals_RHF();
    void sort_OOOV_integrals_RHF();
    void run_simult_dct_oo_RHF();
    void run_simult_dct_RHF();
    void process_so_ints_RHF();
    void build_cumulant_intermediates_RHF();
    void form_density_weighted_fock_RHF();
    void compute_F_intermediate_RHF();
    double compute_cumulant_residual_RHF();
    void update_cumulant_jacobi_RHF();
    void compute_scf_energy_RHF();
    void compute_dct_energy_RHF();
    void compute_orbital_rotation_jacobi_RHF();
    double compute_orbital_residual_RHF();
    void rotate_orbitals_RHF();
    void compute_orbital_gradient_VO_RHF(bool separate_gbargamma);
    void compute_orbital_gradient_OV_RHF(bool separate_gbargamma);
    void compute_unrelaxed_density_VVVV_RHF(bool cumulant_only);
    void compute_unrelaxed_separable_density_VVVV_RHF();
    void compute_unrelaxed_density_OVOV_RHF(bool cumulant_only);
    void compute_unrelaxed_separable_density_OVOV_RHF();
    void compute_unrelaxed_density_OOVV_RHF(bool cumulant_only);
    void compute_unrelaxed_density_OOOO_RHF(bool cumulant_only);
    void compute_unrelaxed_separable_density_OOOO_RHF();
    void compute_gradient_RHF();
    void gradient_init_RHF();
    void compute_lagrangian_OO_RHF(bool separate_gbargamma);
    void compute_lagrangian_VV_RHF(bool separate_gbargamma);
    void compute_ewdm_odc_RHF();
    void print_opdm_RHF();
    void compute_R_AA_and_BB();
    void presort_mo_tpdm_AB();
    void presort_mo_tpdm_AA();
    void construct_oo_density_RHF();

    // UHF-reference DCT
    void compute_gradient_UHF();
    double compute_energy_UHF();

    // Validation
    void validate_energy();
    void validate_opdm();
    void validate_gradient();

    bool augment_b(double* vec, double tol);
    /// Controls convergence of the orbital updates
    bool orbitalsDone_;
    /// Controls convergence of the density cumulant updates
    bool cumulantDone_;
    /// Controls convergence of the DCT energy
    bool energyConverged_;
    /// Whether the user requested the DCT functional that is variationally orbitally-optimized
    bool orbital_optimized_;
    /// Whether the user requested the DCT functional that computes the non-idempotent part of the OPDM exactly from
    /// the density cumulant
    bool exact_tau_;
    /// The amount of information to print
    int print_;
    /// The number of unique pairs of symmetrized atomic orbitals
    int ntriso_;
    /// The number of active alpha electrons
    int nalpha_;
    /// The number of active beta electrons
    int nbeta_;
    /// The number of virtual alpha orbitals
    int navir_;
    /// The number of virtual beta orbitals
    int nbvir_;
    /// The maximum size of the DIIS subspace
    int maxdiis_;
    /// The number of DIIS vectors needed for extrapolation
    int mindiisvecs_;
    /// The maximum number of iterations
    int maxiter_;
    /// The current number of macroiteration for energy or gradient computation
    int iter_;
    // Quadratically-convergent DCT
    /// The total number of independent pairs for the current NR step
    int nidp_;
    /// The number of orbital independent pairs for the current NR step (Alpha spin)
    int orbital_idp_a_;
    /// The number of orbital independent pairs for the current NR step (Beta spin)
    int orbital_idp_b_;
    /// The total number of orbital independent pairs for the current NR step
    int orbital_idp_;
    /// The number of cumulant independent pairs for the current NR step (Alpha-Alpha spin)
    int cumulant_idp_aa_;
    /// The number of cumulant independent pairs for the current NR step (Alpha-Beta spin)
    int cumulant_idp_ab_;
    /// The number of cumulant independent pairs for the current NR step (Beta-Beta spin)
    int cumulant_idp_bb_;
    /// The total number of cumulant independent pairs for the current NR step
    int cumulant_idp_;
    /// The maximum number of IDPs ever possible
    int dim_;
    /// The maximum number of IDPs possible for orbital updates
    int dim_orbitals_;
    /// The maximum number of IDPs possible for cumulant updates
    int dim_cumulant_;
    /// The lookup array that determines which compound indices belong to orbital IDPs and which don't
    std::vector<bool> lookup_orbitals_;
    /// The lookup array that determines which compound indices belong to cumulant IDPs and which don't
    std::vector<bool> lookup_cumulant_;
    /// The number of the guess subspace vectors for the Davidson diagonalization
    int nguess_;
    /// The dimension of the subspace in the Davidson diagonalization
    int b_dim_;
    /// Convergence of the residual in the Davidson diagonalization
    double r_convergence_;
    /// The number of vectors that can be added during the iteration in the Davidson diagonalization
    int n_add_;
    /// The number of eigenvalues requested in the stability check
    int nevals_;
    /// The maximum size of the subspace in stability check
    int max_space_;
    /// The number of occupied alpha orbitals per irrep
    Dimension naoccpi_;
    /// The number of occupied beta orbitals per irrep
    Dimension nboccpi_;
    /// The number of virtual alpha orbitals per irrep
    Dimension navirpi_;
    /// The number of virtual beta orbitals per irrep
    Dimension nbvirpi_;
    /// Alpha occupied MO offset
    std::vector<int> aocc_off_;
    /// Alpha virtual MO offset
    std::vector<int> avir_off_;
    /// Beta occupied MO offset
    std::vector<int> bocc_off_;
    /// Beta virtual MO offset
    std::vector<int> bvir_off_;
    /// The nuclear repulsion energy in Hartree
    double enuc_;
    /// The cutoff below which and integral is assumed to be zero
    double int_tolerance_;
    /// The RMS value of the error vector after the SCF iterations
    double orbitals_convergence_;
    /// The RMS value of the change in lambda after the lambda iterations
    double cumulant_convergence_;
    /// The RMS value of the change in the orbital response
    double orbital_response_rms_;
    /// The RMS value of the change in the cumulant response
    double cumulant_response_rms_;
    /// The RMS value of the change in the coupling of orbital and cumulant response
    double response_coupling_rms_;
    /// The convergence criterion for the lambda iterations
    double cumulant_threshold_;
    /// The convergence criterion for the scf iterations
    double orbitals_threshold_;
    /// The convergence criterion for energy
    double energy_threshold_;
    /// The convergence that must be achieved before DIIS extrapolation starts
    double diis_start_thresh_;
    /// The SCF component of the energy
    double scf_energy_;
    /// The Lambda component of the energy
    double lambda_energy_;
    /// The previous total energy
    double old_total_energy_;
    /// The updated total energy
    double new_total_energy_;
    /// The Tikhonow regularizer used to remove singularities (c.f. Taube and Bartlett, JCP, 2009)
    double regularizer_;
    /// Level shift for denominators in orbital updates
    double orbital_level_shift_;
    /// The threshold for the norm of the residual part of the subspace (|b'> = |b'> - |b><b|b'>) that is used to
    /// augment the subspace
    double vec_add_tol_;
    /// Level shift applied to the diagonal of the density-weighted Fock operator
    double energy_level_shift_;

    /// The alpha occupied eigenvectors, per irrep
    SharedMatrix aocc_c_;
    /// The beta occupied eigenvectors, per irrep
    SharedMatrix bocc_c_;
    /// The alpha virtual eigenvectors, per irrep
    SharedMatrix avir_c_;
    /// The beta virtual eigenvectors, per irrep
    SharedMatrix bvir_c_;
    /// The Tau matrix in the AO basis, stored by irrep, to perturb the alpha Fock matrix
    SharedMatrix tau_so_a_;
    /// The Tau matrix in the AO basis, stored by irrep, to perturb the beta Fock matrix
    SharedMatrix tau_so_b_;
    /// The Tau matrix in the MO basis (alpha occupied)
    Matrix aocc_tau_;
    /// The Tau matrix in the MO basis (beta occupied)
    Matrix bocc_tau_;
    /// The Tau matrix in the MO basis (alpha virtual)
    Matrix avir_tau_;
    /// The Tau matrix in the MO basis (beta virtual)
    Matrix bvir_tau_;
    /// The perturbed Tau matrix in the MO basis (alpha occupied)
    Matrix aocc_ptau_;
    /// The perturbed Tau matrix in the MO basis (beta occupied)
    Matrix bocc_ptau_;
    /// The perturbed Tau matrix in the MO basis (alpha virtual)
    Matrix avir_ptau_;
    /// The perturbed Tau matrix in the MO basis (beta virtual)
    Matrix bvir_ptau_;
    /// The Kappa in the MO basis (alpha occupied)
    SharedMatrix kappa_mo_a_;
    /// The Kappa in the MO basis (beta occupied)
    SharedMatrix kappa_mo_b_;
    /// The overlap matrix in the AO basis
    SharedMatrix ao_s_;
    /// The one-electron integrals in the SO basis
    Matrix so_h_;
    /// The alpha Fock matrix in the SO basis
    SharedMatrix Fa_;
    /// The beta Fock matrix in the SO basis
    SharedMatrix Fb_;
    /// The alpha Fock matrix in the MO basis
    SharedMatrix moFa_;
    /// The beta Fock matrix in the MO basis
    SharedMatrix moFb_;
    /// The alpha density-weighted Fock matrix in the MO basis
    SharedMatrix Ftilde_a_;
    /// The beta density-weighted Fock matrix in the MO basis
    SharedMatrix Ftilde_b_;
    /// The inverse square root overlap matrix in the SO basis
    SharedMatrix s_half_inv_;
    /// The old full alpha MO coefficients
    SharedMatrix old_ca_;
    /// The old full beta MO coefficients
    SharedMatrix old_cb_;
    /// The alpha kappa matrix in the SO basis
    SharedMatrix kappa_so_a_;
    /// The beta kappa matrix in the SO basis
    SharedMatrix kappa_so_b_;
    /// The alpha SCF error vector
    SharedMatrix scf_error_a_;
    /// The beta SCF error vector
    SharedMatrix scf_error_b_;
    // Quadratically-convergent DCT
    /// The orbital gradient in the MO basis (Alpha spin)
    SharedMatrix orbital_gradient_a_;
    /// The orbital gradient in the MO basis (Beta spin)
    SharedMatrix orbital_gradient_b_;
    /// Orbital and cumulant gradient in the basis of IDP
    SharedVector gradient_;
    /// Contribution of the Fock matrix to the diagonal part of the Hessian. Used as preconditioner for conjugate
    /// gradient procedure
    SharedVector Hd_;
    /// The step vector in the IDP basis
    SharedVector X_;
    /// Sigma vector in the basis of IDP (the product of the off-diagonal part of the Hessian with the step vector X)
    SharedVector sigma_;
    /// The conjugate direction vector in the IDP basis for conjugate gradient procedure
    SharedVector D_;
    /// The residual vector in the IDP basis for conjugate gradient procedure
    SharedVector R_;
    /// The search direction vector in the IDP basis for conjugate gradient procedure
    SharedVector S_;
    /// The new element of Krylov subspace vector in the IDP basis for conjugate gradient procedure
    SharedVector Q_;
    /// The subspace vector in the Davidson diagonalization procedure
    SharedMatrix b_;
    /// Orbital parameters. Specifically, the generators of the orbital rotations with respect to
    /// the reference orbitals. Dimension: nmo_ x nmo_.
    SharedMatrix Xtotal_a_;
    SharedMatrix Xtotal_b_;

    /// Used to align things in the output
    std::string indent;

    // Density-Fitting DCT
    /// Set DF variables. Print header.
    void initialize_df();
    /// Calculate memory required for density-fitting
    void estimate_df_memory() const;
    /// Construct b(Q|mn) = Sum_P (mn|P) [J^-1/2]_PQ
    void build_df_b();
    /// Form AO basis b(Q|mu,nu)
    Matrix formb_ao(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                    std::shared_ptr<BasisSet> zero, const Matrix& Jm12, const std::string& name);
    /// Form J(P|Q)^-1/2 and write to disk
    Matrix formJm12(std::shared_ptr<BasisSet> auxiliary, const std::string& name);
    /// Transform SO-basis b(Q, mn) to MO-basis b(Q, pq)
    void transform_b_so2mo();
    /// Transform b(Q|mu,nu) from AO basis to SO basis
    DFTensor transform_b_ao2so(const Matrix& bQmn_ao) const;
    /// Transform b(Q|mu,nu) from AO basis to SO basis
    Matrix transform_b_so2ao(const DFTensor& bQmn_so) const;
    /// Build density-fitted <VV||VV>, <vv||vv>, and <Vv|Vv> tensors in G intermediate
    void build_DF_tensors_RHF();
    void build_DF_tensors_UHF();
    void construct_metric_density(const std::string& basis_type);
    void three_idx_cumulant_density();
    void three_idx_cumulant_density_RHF();
    void three_idx_separable_density();
    DFTensor three_idx_cumulant_helper(DFTensor& temp, const Matrix& J, const Matrix& bt1, const Matrix& bt2);
    DFTensor three_idx_separable_helper(const Matrix& Q, const Matrix& J, const Matrix& RDM, const Matrix& C_subset);

    /// Form density-fitted MO-basis TEI g(OV|OV) in chemists' notation
    void form_df_g_ovov();
    /// Form density-fitted MO-basis TEI g(OO|OO) in chemists' notation
    void form_df_g_oooo();
    /// Form density-fitted MO-basis TEI g(VV|OO) in chemists' notation
    void form_df_g_vvoo();
    /// Form density-fitted MO-basis TEI g(VO|OO) in chemists' notation
    void form_df_g_vooo();
    /// Form density-fitted MO-basis TEI g(OV|VV) in chemists' notation
    void form_df_g_ovvv();
    /// Form density-fitted MO-basis TEI g<VV|VV> in physists' notation
    void form_df_g_vvvv();
    /// Form MO-based Gbar*Gamma
    void build_gbarGamma_RHF();
    void build_gbarGamma_UHF();
    /// Form gbar<ab|cd> * lambda <ij|cd>
    void build_gbarlambda_RHF_v3mem();
    void build_gbarlambda_UHF_v3mem();

    // Density-Fitting DCT
    /// Auxiliary basis
    std::shared_ptr<BasisSet> auxiliary_;
    /// Auxiliary basis for SCF terms in DCT
    std::shared_ptr<BasisSet> auxiliary_scf_;
    /// Primary basis
    std::shared_ptr<BasisSet> primary_;
    /// Number of total primary basis functions
    int nn_;
    /// Number of total auxilliary basis functions
    int nQ_;
    int nQ_scf_;
    /// Number of alpha occupied orbitals
    int naocc_;

    /// b(Q|mu,nu)
    DFTensor bQmn_so_;
    DFTensor bQmn_so_scf_;
    /// b(Q|i, j)
    DFTensor bQijA_mo_;
    DFTensor bQijB_mo_;
    /// b(Q|i, a)
    DFTensor bQiaA_mo_;
    DFTensor bQiaB_mo_;
    /// b(Q|a, b)
    DFTensor bQabA_mo_;
    DFTensor bQabB_mo_;

    /// The Tau in the MO basis (All)
    Matrix mo_tauA_;
    Matrix mo_tauB_;
    /// MO-based (Gbar Tau + Gbar Kappa)
    Matrix mo_gbarGamma_A_;
    Matrix mo_gbarGamma_B_;
    /// MO-based Gamma <r|s>
    Matrix mo_gammaA_;
    Matrix mo_gammaB_;

    std::map<std::string, Slice> slices_;
};

}  // namespace dct
}  // namespace psi

#endif  // Header guard
