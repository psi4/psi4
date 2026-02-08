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

#include <algorithm>
#include <cmath>

#include "psi4/physconst.h"

#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libqt/qt.h"

#include "cghf.h"

#include <Einsums/LinearAlgebra.hpp>
#include <Einsums/TensorAlgebra.hpp>

namespace psi {
namespace scf {

CGHF::CGHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {
    common_init();
}

CGHF::CGHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
           std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio) {
    common_init();
}

CGHF::~CGHF() {}

// Define global einsums::BlockTensors and variables needed throughout the CGHF class
// Nathan: Hey Matt, I was getting weird timer issues because you call things inside here like
// Form_H and find_occupation. This should be done later.
void CGHF::common_init() {
    name_ = "CGHF";

    /*** Begin useless definitions ***/
    // Afaik, these are empty matrices and have no effect on CGHF.
    // However, to avoid a segmentation fault, we must create them.
    Ca_ = SharedMatrix(factory_->create_matrix("alpha MO coefficients (C)"));
    Cb_ = SharedMatrix(factory_->create_matrix("beta MO coefficients (C)"));
    Fa_ = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_ = SharedMatrix(factory_->create_matrix("F beta"));
    Da_ = SharedMatrix(factory_->create_matrix("SCF alpha density"));
    Db_ = SharedMatrix(factory_->create_matrix("SCF beta density"));
    Va_ = SharedMatrix(factory_->create_matrix("G alpha"));
    Vb_ = SharedMatrix(factory_->create_matrix("G beta"));
    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_a_->set_name("alpha orbital energies");
    epsilon_b_ = SharedVector(factory_->create_vector());
    epsilon_b_->set_name("beta orbital energies");
    same_a_b_dens_ = false;
    same_a_b_orbs_ = false;
    /*** End useless definitions ***/

    nuclearrep_ = molecule_->nuclear_repulsion_energy({0.0, 0.0, 0.0});

    // => DIIS variables <=
    err_vecs = std::deque<einsums::BlockTensor<std::complex<double>, 2>>(0);
    Fdiis = std::deque<einsums::BlockTensor<std::complex<double>, 2>>(0);
    diis_coeffs = std::vector<std::complex<double>>(0);
    error_doubles = std::vector<std::complex<double>>(0);

    // Sets nalphapi_ and nbetapi_
    find_occupation();

    form_H();         // Fills a single spin block (AA or BB) with the core Hamiltonian

    // Combine nalphapi_ and nbetapi_ to give nelecpi_ for forming the density matrix
    for (int h = 0; h < nirrep_; h++) {
        nelecpi_.push_back(nalphapi_[h] + nbetapi_[h]);
    };

    // => GLOBAL VARS <=
    // here since Einsums requires a vector
    // G_mat: two-electron ERI as a SharedMatrix
    // nelecpi_: nalphapi_ + nbetapi_
    G_mat = mintshelper()->ao_eri();

    // Double size of nsopi_ for GHF
    for (int h = 0; h < nirrep_; h++) {
        irrep_sizes_.push_back(2 * nsopi_[h]);
    };

    // => BLOCKTENSORS <= will be (2nsopi_ x 2nsopi_) unless listed otherwise
    // Note that all of these are complex except for EINS_, EINX_, and the Fock evals
    // Core Hamiltonian, Fock, Orthogonalized Fock, and Coefficient Matrices
    F0_ = einsums::BlockTensor<std::complex<double>, 2>("F0", irrep_sizes_);
    F_ = einsums::BlockTensor<std::complex<double>, 2>("F", irrep_sizes_);
    Fp_ = einsums::BlockTensor<std::complex<double>, 2>("EINX_", irrep_sizes_);
    C_ = einsums::BlockTensor<std::complex<double>, 2>("C", irrep_sizes_);

    // Gradient FDS - SDF
    FDSmSDF_ = einsums::BlockTensor<std::complex<double>, 2>("C", irrep_sizes_);

    // Spin blocked overlap matrix and orthogonalization matrix.
    // Needed as metric to compute gradient with FDS - SDF
    EINS_ = einsums::BlockTensor<double, 2>("Ovlp", irrep_sizes_);
    EINX_ = einsums::BlockTensor<std::complex<double>, 2>("Orth", irrep_sizes_);

    // Eigenvalues and eigenvectors from diagonalized Fock matrix
    Fevecs_ = einsums::BlockTensor<std::complex<double>, 2>("Fock eigenvectors", irrep_sizes_);
    Fevals_ = einsums::BlockTensor<double, 1>("Fock evals", irrep_sizes_);

    // Density, coulomb, exchange matrices
    D_ = einsums::BlockTensor<std::complex<double>, 2>("D", irrep_sizes_);
    J_ = einsums::BlockTensor<std::complex<double>, 2>("J", irrep_sizes_);
    K_ = einsums::BlockTensor<std::complex<double>, 2>("K", irrep_sizes_);

    // Temporary matrices for intermediate steps
    temp1_ = einsums::BlockTensor<std::complex<double>, 2>("temp1", irrep_sizes_);
    temp2_ = einsums::BlockTensor<std::complex<double>, 2>("temp2", irrep_sizes_);

    F0_.zero();
    F_.zero();
    EINS_.zero();
    EINX_.zero();
    C_.zero();
    Fevecs_.zero();
    Fevals_.zero();
    D_.zero();
    J_.zero();
    K_.zero();
    temp1_.zero();
    temp2_.zero();
    FDSmSDF_.zero();

    subclass_init();  // This appears to set up DFT stuff and external potentials

    preiterations();  // CGHF preiterations() routine
}

// Needed for initializing the Einsums spin-blocked overlap matrix EINS_
// and subsequently the orthogonalization matrix EINX_
// Also builds the spin-blocked Hamiltonian as F0_
void CGHF::preiterations() {
    // => FILL SPIN BLOCKED OVERLAP AND HAMILTONIAN MATRICES <=
    //
    // Note that nsopi_[h] will be HALF of irrep_sizes_[h]
    // EINS_ and F0_ filled within the same loop
    for (int h = 0; h < nirrep_; h++)
        for (int p = 0; p < nsopi_[h]; p++)
            for (int q = 0; q < nsopi_[h]; q++) {
                // Offset by nsopi_[h] to get the beta index for p and q
                int beta_p = p + nsopi_[h];
                int beta_q = q + nsopi_[h];

                EINS_[h].subscript(p, q) = S_->get(h, p, q);  // overlap AA spin block
                EINS_[h].subscript(beta_p, beta_q) =
                    S_->get(h, p, q);  // overlap BB spin block is same as AA since same spatial orbitals

                F0_[h].subscript(p, q) = H_->get(h, p, q);            // core Hamiltonian AA spin block
                F0_[h].subscript(beta_p, beta_q) = H_->get(h, p, q);  // core Hamiltonian BB spin block
            }

    // Constructs orthogonalization matrix X_ using Lowdin symmetric orthogonalization
    //
    // X = S^{-1/2}
    //
    // where eigenvalues less than the threshold 1e-7 are set to 0, and the eigenvectors removed
    HF::form_Shalf();
    for (int h = 0; h < nirrep_; h++)
        for (int j = 0; j < nsopi_[h]; j++)
            for (int k = 0; k < nsopi_[h]; k++) {
                EINX_[h].subscript(j, k) = X_->get(h, j, k);
                EINX_[h].subscript(j + nsopi_[h], k + nsopi_[h]) = X_->get(h, j, k);
            }

    // rotated_sad_guess(); // TODO modify the existing Psi4 SAD scheme (UHF, RHF, etc) since it causes errors otherwise
    //  e.g. the only way for this to work is by using a CORE guess

    if (options_.get_str("GUESS") == "CORE") {
        // form_F();
    }
}

void CGHF::finalize() {
    HF::finalize();
}

void CGHF::form_V() {}

/*
 * Does an explicit 4-index for loop to compute JKwK_
 *
 * The density matrix is decomposed into its 4 spin-blocks and evaluated separately.
 * According to the Slater-Condon rules, the off-diagonal elements (alpha-beta, beta-alpha)
 * have zero Coulomb energy contributions. Therefore, only J_aa and J_bb need to be evaluated
 *
 * There is only one symmetry to take advantage of: K_ab = -K_ba† (negative adjoints)
 * Thus, 5 terms are computed here:
 *
 * J_aa/J_bb =      D_sr, G_pqrs -> J_pq
 *
 * K_aa/K_bb/K_ab = D_sr, G_psrq -> K_pq
 *
 * which are of course computed using their respective density spin block (e.g. K_aa is D_aa)
 *
 */

void CGHF::form_G() {
    J_.zero();
    K_.zero();

    int nso = nso_;
    int dim = 2 * nso;

    // TODO check G_mat dimensions for irreps

    for (int h = 0; h < nirrep_; h++)
        for (int p = 0; p < nsopi_[h]; p++)
            for (int q = 0; q < nsopi_[h]; q++) {
                for (int s = 0; s < nsopi_[h]; s++) {
                    for (int r = 0; r < nsopi_[h]; r++) {
                        auto g_pqrs = G_mat->get(h, p * nsopi_[h] + q, r * nsopi_[h] + s);
                        auto g_psrq = G_mat->get(h, p * nsopi_[h] + s, r * nsopi_[h] + q);

                        auto D_AA = D_[h].subscript(s, r);
                        auto D_BB = D_[h].subscript(s + nsopi_[h], r + nsopi_[h]);

                        auto sum_D = D_AA + D_BB;

                        J_[h].subscript(p, q) += g_pqrs * sum_D - g_psrq * D_AA;
                        J_[h].subscript(p + nsopi_[h], q + nsopi_[h]) += g_pqrs * sum_D - g_psrq * D_BB;
                    }
                }
            }
}

/* F = H + J - K */
void CGHF::form_F() {
    F_.zero();
    F_ += F0_;
    F_ += J_;
    F_ -= K_;
}

/*
 * 1. Orthogonalizes the Fock matrix via Fp_ = X_dagger * F * X_
 * 2. Then diagonalizes Fp_ to give a eigenpairs in the MO basis
 * 3. the eigenvectors are back-transformed back to the AO basis
 */
void CGHF::form_C(double shift) {
    if (shift != 0.0) {
        throw PSIEXCEPTION("CGHF does not support E-shift.");
    }

    // => ORTHOGONALIZE FOCK <=
    // F_ @ EINX_ = temp1_

    if (options_.get_bool("DIIS") && iteration_ > options_.get_int("DIIS_START")) {
        auto error_trace = do_diis();
    } else {
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, EINX_, std::complex<double>{0.0},
                                                    &temp1_);

        // EINX_ @ temp1 = Fp_
        einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp1_, std::complex<double>{0.0},
                                                   &Fp_);
    }

    // => DIAGONALIZE FOCK <=
    // Fevecs_ = Fp_;

    // Must be done in a loop over irreps -- cannot just do it with the BlockTensor
    // Nathan: Why?
    for (int h = 0; h < nirrep_; h++) {
        // Do not diagonalize 0x0 matrix
        if (nsopi_[h] == 0) continue;

        // Hermitian eigensolver
        einsums::linear_algebra::heev<true>(&Fp_[h], &Fevals_[h]);

        // heev retuns the wrong side, so we need to take the inverse (hermitian adjoint)
        // i.e. einsums::linear_algebra::invert(&Fp_[h]);
        // TODO: rework math so we don't have to do this step!
        // Also, this can be done with a blas call

        temp1_[h].zero();
        for (int p = 0; p < irrep_sizes_[h]; p++)
            for (int q = 0; q < irrep_sizes_[h]; q++) {
                auto Fp_pq = Fp_[h].subscript(p, q);

                temp1_[h].subscript(q, p) = std::conj(Fp_pq);
            }
        Fp_[h] = temp1_[h];
    }

    // => BACK TRANSFORM <=
    // EINX_ @ Fevecs_ = C_
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, EINX_, Fp_, std::complex<double>{0.0}, &C_);

    find_occupation();
}

void CGHF::compute_SAD_guess(bool natorb) {
    // Fills the density matrix with the precomputed SAD guess from sad.cc
    // This seems to be more robust than SAP at the moment since SAP
    HF::compute_SAD_guess(natorb);

    for (int h = 0; h < nirrep_; h++)
        for (int p = 0; p < nsopi_[h]; p++)
            for (int q = 0; q < nsopi_[h]; q++) {
                auto Da_val = Da_->get(h, p, q);
                auto Db_val = Db_->get(h, p, q);

                D_[h].subscript(p, q) = Da_val;
                D_[h].subscript(p + nsopi_[h], q + nsopi_[h]) = Db_val;
            }

    temp2_.zero();

    temp1_.zero();

    for (int h = 0; h < nirrep_; h++)
        for (int p = 0; p < irrep_sizes_[h]; p++)
            for (int q = 0; q < irrep_sizes_[h]; q++) {
                temp1_[h].subscript(p, q) = EINS_[h].subscript(p, q);
            }

    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, D_, temp1_, std::complex<double>{0.0},
                                                &temp2_);
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, temp2_, temp2_, std::complex<double>{0.0},
                                               &D_);
}

/*
 * Fills the occupied coefficient matrix Cocc_ as well as the complex conjugate cCocc_
 * Then constructs the density matrix D_ with D_uv = C_ui * C_vi_{conj}
 */
void CGHF::form_D() {
    D_.zero();
    temp1_.zero();
    temp2_.zero();

    // Simply fills Cocc_ and cCocc_ with the (2*nsopi_ x nelecpi_) matrix
    // Note that Cocc_ and cCocc_ are both (2*nsopi_ x 2*nsopi_), but the 'remainder'
    // are all zeros and contribute nothing
    // This is just a minor inefficiency, since these matrices are larger than they should be
    for (int h = 0; h < nirrep_; h++) {
        for (int j = 0; j < irrep_sizes_[h]; j++) {
            for (int k = 0; k < nelecpi_[h]; k++) {
                auto C_jk = C_[h].subscript(j, k);
                temp1_[h].subscript(j, k) = C_jk;
                temp2_[h].subscript(j, k) = std::conj(C_jk);
            }
        }
    }

    // Performs einsums contraction ui,vi->uv with temp1_, temp2_ -> D_)
    einsums::tensor_algebra::einsum(einsums::Indices{einsums::index::u, einsums::index::v}, &D_,     // D_uv
                                    einsums::Indices{einsums::index::u, einsums::index::i}, temp1_,  // Cocc_ui
                                    einsums::Indices{einsums::index::v, einsums::index::i}, temp2_   // Cocc.conj.T_vi
    );
}

void CGHF::damping_update(double damping_percentage) {}

/*
 * E = E_1e + E_2e + E_nuc
 *
 * E_1e = trace(D_ • F0_)
 * E_2e = trace(D_ • 0.5*JK_)
 */
double CGHF::compute_E() {
    double kinetic_E = 0.0;
    double one_electron_E = 0.0;
    double two_E = 0.0;
    double exchange_E = 0.0;

    // Because Psi4 likes these energies (me too) decomposed, it seems both easier and more
    // efficient to just knock them out in the same for loop
    for (int h = 0; h < nirrep_; h++) {
        if (nsopi_[h] == 0) continue;

        for (int i = 0; i < irrep_sizes_[h]; i++)
            for (int j = 0; j < irrep_sizes_[h]; j++) {
                auto Dji = D_[h].subscript(j, i);

                one_electron_E += (F0_[h].subscript(i, j) * Dji).real();  // F0_ij * D_ji
                two_E += (J_[h].subscript(i, j) * Dji).real();            // J_ij * D_ji
            }
    }

    // Add these to the global energies struct
    // These are seen in the output.dat file
    energies_["Nuclear"] = nuclearrep_;
    energies_["Kinetic"] = 1.0;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = 0.5 * two_E;

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += 0.5 * two_E;

    // outfile->Printf("  Core energy : %.15f\n", one_electron_E);
    // outfile->Printf("  JK energy   : %.15f\n", two_E);
    // outfile->Printf("  Total energy: %.15f\n", Etotal);
    return Etotal;
}

/*
 * Simply computes the 1e energy
 */
double CGHF::compute_initial_E() {
    double one_electron_E = 0.0;

    for (int h = 0; h < nirrep_; h++)
        for (int i = 0; i < irrep_sizes_[h]; i++)
            for (int j = 0; j < irrep_sizes_[h]; j++) {
                auto Dji = D_[h].subscript(j, i);
                one_electron_E += (F0_[h].subscript(i, j) * Dji).real();  // F0_ij * D_ji
            }

    // outfile->Printf("  Core energy: %.15f\n", one_electron_E);
    return one_electron_E + nuclearrep_;
}

std::complex<double> CGHF::do_diis() {
    int diis_max = options_.get_int("DIIS_MAX_VECS");
    int diis_count = 0;

    // FDS-SDF
    form_FDSmSDF();

    auto ortho_error =
        einsums::BlockTensor<std::complex<double>, 2>("Orthogonalized FDSmSDF", irrep_sizes_);  // eorth_it
    auto temp1 = einsums::BlockTensor<std::complex<double>, 2>("Temp FDSmSDF", irrep_sizes_);
    auto ecurr = einsums::BlockTensor<std::complex<double>, 2>("Current error", irrep_sizes_);

    temp1_.zero();
    ortho_error.zero();
    ecurr.zero();
    temp2_.zero();

    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, EINX_, std::complex<double>{0.0},
                                                &temp2_);
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp2_, std::complex<double>{0.0},
                                               &Fp_);

    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, FDSmSDF_, EINX_, std::complex<double>{0.0},
                                                &temp1_);
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp1_, std::complex<double>{0.0},
                                               &ortho_error);

    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, ortho_error, ortho_error,
                                               std::complex<double>{0.0}, &ecurr);

    auto error_trace = std::complex<double>{0.0, 0.0};

    // Get collective trace of all irreps
    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < irrep_sizes_[i]; j++) {
            error_trace = ecurr[i](j, j);
        }
    }

    auto abs_trace = std::abs(error_trace);

    // If we have enough error vectors in storage, check for the worst one to remove
    if (err_vecs.size() == diis_max) {
        // Find the one with the max error to eliminate from storage
        // Give it some ludicrously low error so that anything is larger
        // C++ isn't like Python where we can assign the var in the loop
        // at least not at first
        auto max_error = std::complex<double>{-10000000000, -10000000000};
        int max_error_ind = 0;

        for (int i = 0; i < diis_max; i++) {
            auto curr_error = error_doubles.at(i);

            // You can't compare two complex numbers directly, so I compare the real parts
            // Doing it with imaginary is very tricky since these are very very close to 0
            // When I checked imag as well, we were converging at 17 iterations instead of 14
            if (curr_error.real() > max_error.real()) {
                max_error = curr_error;
                max_error_ind = i;
            }
        }
        /*
         * Use this to replace the Fock matrix with the most error
         * Seems to be more applicable to what I'm doing
         */
        // Then make the replacement based off the index
        Fdiis.at(max_error_ind) = Fp_;
        err_vecs.at(max_error_ind) = FDSmSDF_;

        // Norm for the replacement error vector, which is why the index is different here compared to below
        auto norm = einsums::linear_algebra::dot(err_vecs.at(max_error_ind), err_vecs.at(max_error_ind));
        error_doubles.at(max_error_ind) = norm;  // This is set to real for now, not sure what to do about this

        /*
         * Use this to reset the subspace at every nth iteration
         * Might be more useful for generally converging
         * energies rather than longer geometries
         *
        Fdiis.clear();
        err_vecs.clear();
        Fdiis.push_back(Fp_);
        err_vecs.push_back(error_);
        */
    }

    // If not, add the Fock matrix and error matrix to memory
    else {
        err_vecs.push_back(FDSmSDF_);
        Fdiis.push_back(Fp_);
        // Calculate the norm for the last error vector
        auto norm = einsums::linear_algebra::dot(err_vecs.at(err_vecs.size() - 1), err_vecs.at(err_vecs.size() - 1));
        error_doubles.push_back(norm);  // See above comment with norm
    }

    // Next form the 'base' B matrix
    auto B_ = einsums::Tensor<std::complex<double>, 2>("Overlap matrix", err_vecs.size() + 1, err_vecs.size() + 1);
    B_.zero();

    // I do it in a weird way by filling the bottom row and
    // the right column with the 1.0+0i
    for (int i = 0; i < err_vecs.size() + 1; i++) {
        B_(i, err_vecs.size()) = std::complex<double>{1.0, 0.0};
        B_(err_vecs.size(), i) = std::complex<double>{1.0, 0.0};
    }

    // Then simply zero out the last element in the matrix
    B_(err_vecs.size(), err_vecs.size()) = std::complex<double>{0.0, 0.0};

    // Actually populate it with the overlaps

    for (int i = 0; i < err_vecs.size(); i++) {
        for (int j = 0; j < err_vecs.size(); j++) {
            B_(i, j) = {0.0, 0.0};
            for (int p = 0; p < irrep_sizes_[0]; p++) {
                for (int q = 0; q < irrep_sizes_[0]; q++) {
                    B_(i, j) += std::conj(err_vecs[i](q, p)) * err_vecs[j](p, q);
                }
            }
        }
    }

    // Container for the coefficients after gesv
    // This is a weird artifact of gesv. This should be a column vector
    // But it doesn't accept the type? IDK
    auto C_temp = einsums::Tensor<std::complex<double>, 2>("Temp coefficient matrix", 1, err_vecs.size() + 1);

    C_temp.zero();

    for (int i = 0; i < err_vecs.size() + 1; i++) {
        C_temp(0, i) = std::complex<double>{1.0, 0.0};
    }

    einsums::linear_algebra::gesv(&B_, &C_temp);

    // Increase the size of the subspace
    diis_coeffs.resize(err_vecs.size());

    for (int i = 0; i < err_vecs.size(); i++) {
        diis_coeffs.at(i) = C_temp(0, i);
    }

    // Then we can extrapolate the Fock matrix!
    Fp_.zero();

    if (err_vecs.size() == diis_max) {
        for (int i = 0; i < diis_coeffs.size(); i++) {
            einsums::linear_algebra::axpy(diis_coeffs[i], Fdiis[i], &Fp_);
        }
    } else {
        temp2_.zero();
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, F_, EINX_, std::complex<double>{0.0},
                                                    &temp2_);
        einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, EINX_, temp2_, std::complex<double>{0.0},
                                                   &Fp_);
    }

    return error_trace;
}

void CGHF::form_FDSmSDF() {
    auto SComplex = einsums::BlockTensor<std::complex<double>, 2>("Complex overlap matrix", irrep_sizes_);
    SComplex.zero();

    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < irrep_sizes_[i]; j++) {
            for (int k = 0; k < irrep_sizes_[i]; k++) {
                SComplex[i].subscript(j, k) = EINS_[i].subscript(j, k);
            }
        }
    }

    for (int i = 0; i < nirrep_; i++) {
        for (int j = 0; j < irrep_sizes_[i]; j++) {
            for (int k = 0; k < irrep_sizes_[i]; k++) {
                FDSmSDF_[i](j, k) = 0.;
                for (int p = 0; p < irrep_sizes_[i]; p++) {
                    for (int q = 0; q < irrep_sizes_[i]; q++) {
                        FDSmSDF_[i](j, k) += F_[i](j, p) * D_[i](p, q) * SComplex[i](q, k);
                        FDSmSDF_[i](j, k) -= std::conj(F_[i](k, p) * D_[i](p, q) * SComplex[i](q, j));
                    }
                }
            }
        }
    }
}

double CGHF::compute_Dnorm() {
    double dnorm = 0.0;

    for (int h = 0; h < nirrep_; h++) {
        dnorm += einsums::linear_algebra::norm(einsums::linear_algebra::Norm::Frobenius, FDSmSDF_[h]);
    }

    return dnorm;
}

void CGHF::setup_potential() {}

void CGHF::openorbital_scf() {}

std::shared_ptr<CGHF> CGHF::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    auto wfn = Wavefunction::c1_deep_copy(basis);
    // auto hf_wfn = std::make_shared<CGHF>(wfn, functional_, wfn->options(), wfn->psio());

    auto hf_wfn = std::shared_ptr<CGHF>(new CGHF(wfn, functional_, wfn->options(), wfn->psio()));

    return hf_wfn;
}

}  // namespace scf
}  // namespace psi
