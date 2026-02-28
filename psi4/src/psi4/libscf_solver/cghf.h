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

#ifndef PSI4_LIBSCF_SOLVER_CGHF
#define PSI4_LIBSCF_SOLVER_CGHF

#include "psi4/libpsio/psio.hpp"
#include "psi4/libfock/v.h"
#include "hf.h"

#ifdef USING_Einsums
#include <Einsums/Config.hpp>
#include <Einsums/Tensor.hpp>

using SharedBlockTensor = std::shared_ptr<einsums::BlockTensor<std::complex<double>, 2>>;
#endif

namespace psi {
namespace scf {

#ifndef USING_Einsums
class CGHF : public HF {
   public:
    CGHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
      :HF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {
        throw PSIEXCEPTION("Psi4 not built with Einsums. CGHF not available.");
    };

    CGHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options, std::shared_ptr<PSIO> psio)
      :HF(ref_wfn, func, options, psio) {
        throw PSIEXCEPTION("Psi4 not built with Einsums. CGHF not available.");
    }

    // Required for export_wavefunction to build
    double compute_Dnorm() {return 0.0;}
    void preiterations() {};
    void form_FDSmSDF() {};

    // Required to build
    std::shared_ptr<UV> potential_;
    std::shared_ptr<VBase> V_potential() const override { return potential_; };

    ~CGHF() override {};
};
#else

class CGHF : public HF {
   public:
    CGHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional);
    CGHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> functional, Options& options,
         std::shared_ptr<PSIO> psio);
    ~CGHF() override;

    // Declares all BlockTensors, and some useful misc variables needed
    void common_init();

    // Constructs the spin-blocked overlap (EINS_), core Hamiltonian (F0_), and orthogonalization (EINX_) matrices
    void preiterations();
    void save_density_and_energy() {}

    // Allow SAP and SAD initial guesses
    void sap_guess();
    void compute_SAD_guess(bool natorb) override;

    // Form orbital gradient FDSmSDF_ = [F, D]
    void form_FDSmSDF();

    // Compute the norm from the orbital gradient as a second test of convergence
    double compute_Dnorm();

    void form_Shalf() override;

    // Empty function for now -- needed for DFT later
    void form_V() override;

    // Computes J and K either explicitly (4-index) or with RI (3-index)
    void form_G() override;

    // Forms Fock matrix F = F0_ + J - K
    void form_F() override;

    // Orthogonalizes then diagonalizes the Fock matrix to form the coefficient matrix C_
    void form_C(double shift = 0.0) override;

    // Form coefficient matrix from Fock guess (SAP, CORE, etc.)
    void form_initial_C() override;

    // Constructs 1-particle density matrix using the occupied coefficients Cocc_
    // and its conjugate. Stored in temp1_.
    void form_D() override;

    // Empty function for now, but for UHF and RHF, scales the density matrix
    void damping_update(double damping_percentage) override;

    // Compute the energy based purely off F0_ = T + V, with no J and K
    double compute_initial_E() override;

    // Compute 1e and 2e energy separately, then combine with nuclear repulsion
    // energy nuclearrep_ to return a total energy.
    double compute_E() override;

    void setup_potential() override {};

    std::shared_ptr<CGHF> c1_deep_copy(std::shared_ptr<BasisSet> basis);

    // Unsure of what these are, but needed otherwise seg fault
    bool same_a_b_orbs() const { return false; }
    bool same_a_b_dens() const { return false; }

    // Empty functions for now (TODO later)
    std::shared_ptr<UV> potential_;
    std::shared_ptr<VBase> V_potential() const override { return potential_; };

    // If DIIS is enabled, this will update the orthogonalized Fock matrix Fp_
    std::complex<double> do_diis();

   protected:
    SharedMatrix V_mat;
    SharedMatrix S_mat;
    SharedMatrix T_mat;
    SharedMatrix G_mat;
    SharedMatrix F_mat;

    // DIIS variables
    // All 4 of these containers have a MAX size of DIIS_MAX_VECS

    // Holds the grabbed Fock matrices to extrapolate
    std::deque<einsums::BlockTensor<std::complex<double>, 2>> Fdiis;
    // Holds FDSmSDF_ at each iteration (orbital gradients)
    std::deque<einsums::BlockTensor<std::complex<double>, 2>> err_vecs;
    std::vector<std::complex<double>> diis_coeffs;    // Holds the coefficients for each Fock matrix in Fdiis
    std::vector<std::complex<double>> error_doubles;  // RMS errors (real)

    double nuclearrep_;  // Nuclear repulsion energy

    // Core Hamiltonian, Fock, Orthogonalized Fock, and Coefficient Matrices
    SharedBlockTensor F0_;       // Core Hamilton F0 = T + V
    SharedBlockTensor F_;        // Non-orthogonal Fockl: F = T + V + J - K
    SharedBlockTensor FDSmSDF_;  // Fock gradient [F,D]
    SharedBlockTensor Fp_;       // Orthogonalized Fock matrix
    SharedBlockTensor C_;        // Coefficient matrix built after back-trasnformation C' = XC

    // NOTE: EINS_ and EINX_ are spin-blocked variants of S_ and X_ from HF, respectively
    SharedBlockTensor EINS_;  // Spin-blocked overlap matrix -- it is never complex
    SharedBlockTensor EINX_;  // Spin-blocked orthogonalization matrix

    SharedBlockTensor Fevecs_;                // Eigenvectors of Fock matrix
    einsums::BlockTensor<double, 1> Fevals_;  // Eigenvalues of Fock matrix

    SharedBlockTensor D_;   // 1-particle density matrix
    SharedBlockTensor JK_;  // Combined Coulomb and exchange matrix

    // ortho_error and ecurr are specific to DIIS
    SharedBlockTensor ortho_error;  // Orthogonalized gradient error
    SharedBlockTensor ecurr;        // Error at current iteration

    // temp1_ and temp2_ are temporary storage containers for intermediate steps
    // Cocc_ and cCocc_ are not preferred as variable names since there doesn't
    // appear to be a reason to permanently store these (temp1_ and temp2_ are used instead)
    SharedBlockTensor temp1_;
    SharedBlockTensor temp2_;

    // Number of spin orbitals per irrep. Cannot use nsopi_ because Einsums requires a vector.
    std::vector<int> irrep_sizes_;  // Since GHF is spin-blocked, each irrep (h) size will be 2*nsopi_[h]
    std::vector<int> nelecpi_;      // Number of electrons per irrep
};
#endif

}  // namespace scf
}  // namespace psi

#endif
