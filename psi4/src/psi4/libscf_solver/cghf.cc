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

#ifdef USING_Einsums
#include <Einsums/LinearAlgebra.hpp>
#include <Einsums/TensorAlgebra.hpp>
#include <Einsums/Runtime.hpp>

using ComplexMatrix = einsums::BlockTensor<std::complex<double>, 2>;
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace scf {

#ifdef USING_Einsums

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
void CGHF::common_init() {
    name_ = "CGHF";

    // ao_eri lacks irreps and we have no JK object yet so we only support C1
    if (nirrep_ > 1) throw PSIEXCEPTION("USE C1 SYMMETRY!");

    // Some of these are used to pass matrices from HF::guess.
    // However, we get segmentation faults when they are not all defined.
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

    // => DIIS variables <=
    err_vecs = std::deque<std::shared_ptr<ComplexMatrix>>(0);
    Fdiis = std::deque<std::shared_ptr<ComplexMatrix>>(0);
    error_doubles = std::vector<double>(0);

    // nelecpi_ is needed to form_C which is called after every guess. Preiterations
    // is called after guess so this needs to go here.
    find_occupation();
    // Combine nalphapi_ and nbetapi_ to give nelecpi_ for forming the density matrix
    for (int h = 0; h < nirrep_; h++) {
        nelecpi_.push_back(nalphapi_[h] + nbetapi_[h]);
    }

    // Double size of nsopi_ for GHF
    for (int h = 0; h < nirrep_; h++) {
        irrep_sizes_.push_back(2 * nsopi_[h]);
    }

    // => BLOCKTENSORS <= will be (2nsopi_ x 2nsopi_) unless listed otherwise
    // Note that all of these are complex except for the Fock evals
    // Core Hamiltonian, Fock, Orthogonalized Fock, and Coefficient Matrices
    F0_ = std::make_shared<ComplexMatrix>("F0", irrep_sizes_);
    F_ = std::make_shared<ComplexMatrix>("F", irrep_sizes_);
    Fp_ = std::make_shared<ComplexMatrix>("Forth", irrep_sizes_);
    C_ = std::make_shared<ComplexMatrix>("C", irrep_sizes_);

    // Gradient FDS - SDF
    FDSmSDF_ = std::make_shared<ComplexMatrix>("DIIS error", irrep_sizes_);

    // Spin blocked overlap matrix and orthogonalization matrix.
    // Needed as metric to compute gradient with FDS - SDF
    EINS_ = std::make_shared<ComplexMatrix>("Ovlp", irrep_sizes_);
    EINX_ = std::make_shared<ComplexMatrix>("Orth", irrep_sizes_);

    // Eigenvalues and eigenvectors from diagonalized Fock matrix
    Fevecs_ = std::make_shared<ComplexMatrix>("Fock eigenvectors", irrep_sizes_);
    Fevals_ = einsums::BlockTensor<double, 1>("Fock evals", irrep_sizes_);

    // Density and combined Coulomb/ exchange matrices
    D_ = std::make_shared<ComplexMatrix>("D", irrep_sizes_);
    JK_ = std::make_shared<ComplexMatrix>("J", irrep_sizes_);

    // Temporary matrices for intermediate steps
    temp1_ = std::make_shared<ComplexMatrix>("temp1", irrep_sizes_);
    temp2_ = std::make_shared<ComplexMatrix>("temp2", irrep_sizes_);

    F0_->zero();
    F_->zero();
    EINS_->zero();
    EINX_->zero();
    C_->zero();
    Fevecs_->zero();
    Fevals_.zero();
    D_->zero();
    JK_->zero();
    temp1_->zero();
    temp2_->zero();
    FDSmSDF_->zero();

    subclass_init();

    // Initialize einsums and turn off logging to stdout.
    const char* ein_argv[4] = {"psi4\0", "--einsums:no-profiler-report\0", "--einsums:log-level\0", "3\0"};
    einsums::initialize(4, ein_argv);
}

/* Needed for initializing the Einsums spin-blocked overlap matrix EINS_ and
 * subsequently the orthogonalization matrix EINX_ Also builds the spin-blocked
 * Hamiltonian as F0_ */
void CGHF::preiterations() {
    nuclearrep_ = molecule_->nuclear_repulsion_energy({0.0, 0.0, 0.0});
    G_mat = mintshelper()->ao_eri();

    size_t nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    // => FILL SPIN BLOCKED OVERLAP AND HAMILTONIAN MATRICES <=
    // Note that nsopi_[h] will be HALF of irrep_sizes_[h]
    // EINS_ and F0_ filled within the same loop
    for (int h = 0; h < nirrep_; h++) {
        int dim = nsopi_[h];
#pragma omp parallel for schedule(auto) num_threads(nthreads)
        for (int p = 0; p < dim; p++) {
            for (int q = 0; q < dim; q++) {
                // Offset by nsopi_[h] to get the beta index for p and q
                int beta_p = p + dim;
                int beta_q = q + dim;

                // EINS_: [[S_, 0],[0, S_]]
                EINS_->block(h)(p, q) = S_->get(h, p, q);
                EINS_->block(h)(beta_p, beta_q) = S_->get(h, p, q);

                // F0_: [[H_, 0],[0, H_]]
                F0_->block(h)(p, q) = H_->get(h, p, q);
                F0_->block(h)(beta_p, beta_q) = H_->get(h, p, q);

                // F_: [[Fa_, 0],[0, Fb_]]
                F_->block(h)(p, q) = Fa_->get(h, p, q);
                F_->block(h)(beta_p, beta_q) = Fb_->get(h, p, q);
            }
        }
    }

    // If SOC, then add Hx, Hy, Hz to F0_
    if (options_.get_bool("SPIN_ORBIT_COUPLING")) {
        outfile->Printf("\n  CGHF: Spin-orbit coupling selected.\n");
        std::vector<SharedMatrix> H_SO(3);
        if (options_.get_str("RELATIVISTIC") == "ZORA") {
            H_SO = mintshelper()->ao_zora_spin_orbit();
        } else {
            throw PSIEXCEPTION("RELATIVISTIC option must be in [\"ZORA\"] when SPIN_ORBIT_COUPLING TRUE");
        }
        std::complex<double> i(0, 1);

        for (int h = 0; h < nirrep_; h++) {
            int nso = nsopi_[h];
            for (int p = 0; p < nso; p++) {
                for (int q = 0; q < nso; q++) {
                    // Note that I multiply by i to get H_SO Hermitian
                    F0_->block(h)(p, q)         += H_SO[2]->get(p,q)*i; // +Hz [0,0]
                    F0_->block(h)(p+nso, q+nso) -= H_SO[2]->get(p,q)*i; // -Hz [1,1]
                    F0_->block(h)(p, q+nso)     += H_SO[0]->get(p,q)*i; // +Hx [0,1]
                    F0_->block(h)(p+nso, q)     += H_SO[0]->get(p,q)*i; // +Hx [1,0]
                    F0_->block(h)(p, q+nso)     += H_SO[1]->get(p,q);   //+iHy [0,1]
                    F0_->block(h)(p+nso, q)     -= H_SO[1]->get(p,q);   //-iHy [1,0]
                }
            }
        }
    }
}

/*
 * Constructs orthogonalization matrix X_ using Lowdin symmetric orthogonalization
 * X = S^{-1/2}
 * Then places the result in an Einsums Matrix EINX_
 */
void CGHF::form_Shalf() {
    HF::form_Shalf();

    size_t nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    for (int h = 0; h < nirrep_; h++) {
        int dim = nsopi_[h];
#pragma omp parallel for schedule(auto) num_threads(nthreads)
        for (int j = 0; j < dim; j++)
            for (int k = 0; k < dim; k++) {
                EINX_->block(h)(j, k) = X_->get(h, j, k);              // AA spin block
                EINX_->block(h)(j + dim, k + dim) = X_->get(h, j, k);  // BB spin block
            }
    }
}

void CGHF::form_V() {}

/*
 * Does an explicit 4-index for loop to compute JKwK_
 *
 * The density matrix is decomposed into its 4 spin-blocks and evaluated separately.
 * According to the Slater-Condon rules, the off-diagonal elements (alpha-beta, beta-alpha)
 * have->zero Coulomb energy contributions. Therefore, only J_aa and J_bb need to be evaluated
 *
 * There is only one symmetry to take advantage of: K_ab = -K_ba† (negative adjoints)
 * Thus, 5 terms are computed here:
 *
 * J_aa/J_bb =      D_sr, G_pqrs -> J_pq
 *
 * K_aa/K_bb/K_ab = D_sr, G_psrq -> K_pq
 *
 * which are of course computed using their respective density spin block (e.g. K_aa is D_aa)
 */
void CGHF::form_G() {
    JK_->zero();

    // TODO check G_mat dimensions for irreps

    size_t nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    for (int h = 0; h < nirrep_; h++) {
        int dim = nsopi_[h];

        auto alpha = einsums::Range{0, dim};
        auto beta  = einsums::Range{dim, 2*dim};

        einsums::TensorView D_AA = D_->block(h)(alpha, alpha);
        einsums::TensorView D_AB = D_->block(h)(alpha, beta);
        einsums::TensorView D_BA = D_->block(h)(beta, alpha);
        einsums::TensorView D_BB = D_->block(h)(beta, beta);

        auto sum_D = D_AA + D_BB;
        auto &JK = JK_->block(h);

        // collapse(2) gives each thread its own (p,q) pair meaning no race conditions on JK_.
#pragma omp parallel for collapse(2) schedule(auto) num_threads(nthreads)
        for (int p = 0; p < dim; p++)
            for (int q = 0; q < dim; q++)
                for (int s = 0; s < dim; s++)
                    for (int r = 0; r < dim; r++) {
                        auto g_pqrs = G_mat->get(h, p * dim + q, r * dim + s);
                        auto g_psrq = G_mat->get(h, p * dim + s, r * dim + q);

                        // Each line will be Coulomb (J) - exchange (K) contributions
                        // NOTE: the off-diagonal AB/BA blocks have no Coulomb contributions, only exchange
                        JK(p, q) += g_pqrs * sum_D(s, r) - g_psrq * D_AA(s, r);
                        JK(p + dim, q + dim) += g_pqrs * sum_D(s, r) - g_psrq * D_BB(s, r);  // BB
                        JK(p, q + dim) -= D_AB(s, r) * g_psrq;  // AB (exchange only)
                        JK(p + dim, q) -= D_BA(s, r) * g_psrq;  // BA (exchange only)
                    }
    }
}

/* F = H + J - K */
void CGHF::form_F() {
    F_->zero();
    (*F_) += (*F0_);
    (*F_) += (*JK_);
}

/*
 * 1. Orthogonalizes the Fock matrix via Fp_ = X_dagger * F * X_
 * 2. Then diagonalizes Fp_ to give a eigenpairs in the MO basis
 * 3. the eigenvectors are back-transformed back to the AO basis
 */
void CGHF::form_C(double shift) {
    // TODO: allow shift as we should be able to do it within one of the gemm calls below by
    // replacing std::complex<double>{0.0} with the shift
    if (shift != 0.0) throw PSIEXCEPTION("CGHF does not support energy shifting.");

    // => ORTHOGONALIZE FOCK <=
    // Fp_ = X_.conj().T @ F_ @ X_

    temp2_->zero();
    // Orthogonalize Fock matrix
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, *F_, *EINX_, std::complex<double>{0.0},
                                                temp2_.get());
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, *EINX_, *temp2_, std::complex<double>{0.0},
                                               Fp_.get());

    if (options_.get_bool("DIIS") && iteration_ > options_.get_int("DIIS_START")) {
        do_diis();
    }

    // => DIAGONALIZE FOCK <=
    // Must be done in a loop over irreps -- cannot just do it with the BlockTensor

    temp1_->zero();
    for (int h = 0; h < nirrep_; h++) {
        // Do not diagonalize 0x0 matrix
        if (nsopi_[h] == 0) continue;

        // Hermitian eigensolver
        einsums::linear_algebra::heev<true>(&Fp_->block(h), &Fevals_[h]);

        for (int i = 0; i < nsopi_[h]; i++) {
            // Putting half the orbitals in alpha half in beta
            epsilon_a_->set(h, i, Fevals_[h].subscript(2 * i));
            epsilon_b_->set(h, i, Fevals_[h].subscript(2 * i + 1));
        }

        // heev retuns the wrong side, so we need to take the inverse (hermitian adjoint)

        // Takes the conjugate transpose of Fp_[h] (e.g. ij -> ji) to give us the proper eigenvectors
        // NOTE: the template parameters <true> states to take the conjugate
        einsums::tensor_algebra::permute<true>(
            std::complex<double>{0.0}, einsums::Indices{einsums::index::i, einsums::index::j}, &temp1_->block(h),
            std::complex<double>{1.0}, einsums::Indices{einsums::index::j, einsums::index::i}, Fp_->block(h));
    }

    // => BACK TRANSFORM <=
    // EINX_ @ Fevecs_ = C_
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, *EINX_, *temp1_, std::complex<double>{0.0},
                                                C_.get());

    find_occupation();
}

/*
 * Certain guesses in HF:guess() put the initial Fock matrix in SharedMatrix Fa_/Fb_.
 * This puts it in BlockTensor F_ then calls form_C() (which expects F_).
 */
void CGHF::form_initial_C() {
    // find_occupation();
    if (!sad_) {
        for (int h = 0; h < nirrep_; h++) {
            for (int p = 0; p < nsopi_[h]; p++) {
                for (int q = 0; q < nsopi_[h]; q++) {
                    F_->block(h)(p, q) = Fa_->get(h, p, q);
                    F_->block(h)(p + nsopi_[h], q + nsopi_[h]) = Fb_->get(h, p, q);
                }
            }
        }
    }

    form_C();
}

/*
 * Fills the density matrix with the precomputed SAD guess from sad.cc
 * This seems to be more robust than SAP at the moment since SAP
 */
void CGHF::compute_SAD_guess(bool natorb) {
    HF::compute_SAD_guess(natorb);

    for (int h = 0; h < nirrep_; h++)
        for (int p = 0; p < nsopi_[h]; p++)
            for (int q = 0; q < nsopi_[h]; q++) {
                D_->block(h)(p, q) = Da_->get(h, p, q);
                D_->block(h)(p + nsopi_[h], q + nsopi_[h]) = Db_->get(h, p, q);
                C_->block(h)(p, 2*q) = Ca_->get(h, p, q);
                C_->block(h)(p + nsopi_[h], 2*q + 1) = Cb_->get(h, p, q);
            }
}

/*
 * Fills the occupied coefficient matrix Cocc_ as well as the complex conjugate cCocc_
 * Then constructs the density matrix D_ with D_uv = C_ui * C_vi_{conj}
 */
void CGHF::form_D() {
    D_->zero();
    temp1_->zero();
    temp2_->zero();


    // Fills temp1_ and temp2_ with the occupied (2*nsopi_ x nelecpi_) and conjugate
    // occupied matrices (e.g. Cocc)
    // Note that temp1_ and temp2_ are both (2*nsopi_ x 2*nsopi_), but the 'remainder'
    // are all->zeros and contribute nothing
    // This is just a minor inefficiency, since these matrices are larger than they should be
    // NOTE: this is done like this because BlockTensors cannot handle rectangular matrices
    for (int h = 0; h < nirrep_; h++) {
        for (int j = 0; j < irrep_sizes_[h]; j++) {
            for (int k = 0; k < nelecpi_[h]; k++) {
                auto C_jk = C_->block(h)(j, k);
                temp1_->block(h)(j, k) = C_jk;
                temp2_->block(h)(j, k) = std::conj(C_jk);
            }
        }
    }

    // D_ = einsums("ui,vi->uv", temp1_, temp2_)
    einsums::tensor_algebra::einsum(einsums::Indices{einsums::index::u, einsums::index::v}, D_.get(),  // D_uv
                                    einsums::Indices{einsums::index::u, einsums::index::i}, *temp1_,   // Cocc_ui
                                    einsums::Indices{einsums::index::v, einsums::index::i}, *temp2_    // Cocc_vi.conj().T
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

    for (int h = 0; h < nirrep_; h++) {
        if (nsopi_[h] == 0) continue;

        for (int i = 0; i < irrep_sizes_[h]; i++)
            for (int j = 0; j < irrep_sizes_[h]; j++) {
                auto Dji = D_->block(h)(j, i);

                one_electron_E += (F0_->block(h)(i, j) * Dji).real();  // F0_ij * D_ji
                two_E += (JK_->block(h)(i, j) * Dji).real();           // J_ij * D_ji
            }


        auto one_electron_E = einsums::Tensor<std::complex<double>, 0>("One electron E");
        auto two_electron_E = einsums::Tensor<std::complex<double>, 0>("Two electron E");

        //einsums::tensor_algebra::einsum(einsums::Indices{}, &one_electron_E,  // D_uv
        //                                einsums::Indices{einsums::index::i, einsums::index::j}, F0_->block(h),   // Cocc_ui
        //                                einsums::Indices{einsums::index::j, einsums::index::i}, D_->block(h)    // Cocc_vi.conj()
        //);

        //net_one_electron_E += one_electron_E.real();


    }

    energies_["Nuclear"] = nuclearrep_;
    energies_["Kinetic"] = 1.0;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = 0.5 * two_E;

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += 0.5 * two_E;

    return Etotal;
}

double CGHF::compute_initial_E() {
    // TODO: F0_ is not yet populated
    double one_electron_E = 0.0;

    for (int h = 0; h < nirrep_; h++)
        for (int i = 0; i < irrep_sizes_[h]; i++)
            for (int j = 0; j < irrep_sizes_[h]; j++) {
                auto Dji = D_->block(h)(j, i);
                one_electron_E += (F0_->block(h)(i, j) * Dji).real();  // F0_ij * D_ji
            }

    return one_electron_E + nuclearrep_;
}

/*
 * Direct inversion of the iterative subspace (DIIS) for improved convergence
 *
 * The approach taken here is to accumulate Fock matrices until DIIS_MAX_VECS is reached.
 * Like usual, a new Fock matrix is extrapolated at this point. To obtain new vectors and
 * replace an old one, the Fock matrix with the largest error according to [F, D] is
 * substituted for the Fock matrix at the current iteration
 *
 * NOTE: ADIIS is not yet available
 */
void CGHF::do_diis() {
    int diis_max = options_.get_int("DIIS_MAX_VECS");

    // FDS-SDF
    form_FDSmSDF();

    // e.e for e = [F,D]
    double error = 0;
    for (int h = 0; h < nirrep_; h++) {
        auto o = FDSmSDF_->block(h);
        // Avoid using true_dot until einsums fixes segfault on M-intel runners.
        // error += einsums::linear_algebra::true_dot(o, o).real();
        for (int p = 0; p < irrep_sizes_[h]; p++) {
            for (int q = 0; q < irrep_sizes_[h]; q++) {
                error += (std::conj(o.subscript(p, q)) * o.subscript(p, q)).real();
            }
        }
    }

    // Only extrapolate Fp_ when enough vectors are stored
    if (err_vecs.size() < diis_max) {
        error_doubles.push_back(error);
        Fdiis.push_back(Fp_);
        err_vecs.push_back(FDSmSDF_);

        // We create the tensor here instead of diis_max tensors during common_init
        // to account for future DIIS schemes where e.g. vectors are cleared when
        // condition number is too high and so err_vecs.size() is a meaningful quantity.
        Fp_ = std::make_shared<ComplexMatrix>(*Fdiis.back());
        FDSmSDF_ = std::make_shared<ComplexMatrix>(*err_vecs.back());

        return;
    }

    // Find the index with the largest error to eliminate from storage
    double max_error = 0;
    int max_error_ind = 0;

    for (int i = 0; i < diis_max; i++) {
        double curr_error = error_doubles.at(i);

        if (curr_error > max_error) {
            max_error = curr_error;
            max_error_ind = i;
        }
    }

    // Replace max_error with current error
    error_doubles.at(max_error_ind) = error;
    std::swap(Fdiis[max_error_ind], Fp_);
    std::swap(err_vecs[max_error_ind], FDSmSDF_);

    /* Alternatively, we can do this to just replace the last vector */
    // error_doubles.erase(error_doubles.begin());
    // std::swap(Fdiis[0], Fp_);
    // std::swap(err_vecs[0], FDSmSDF_);
    //
    // error_doubles.push_back(error);

    auto B_ = einsums::Tensor<std::complex<double>, 2>("DIIS Error matrix", diis_max + 1, diis_max + 1);
    B_.zero();

    /* B_ has the form (if diis_max == 3):
     *   /  e  e  e -1  \
     *   |  e  e  e -1  |
     *   |  e  e  e -1  |
     *   \ -1 -1 -1  0  /
     */
    for (int i = 0; i < diis_max; i++) {
        B_(i, diis_max) = -1;
        B_(diis_max, i) = -1;
    }

    // Actually populate it with the errors
    for (int i = 0; i < diis_max; i++) {
        for (int j = i; j < diis_max; j++) {
            auto ei = err_vecs[i];
            auto ej = err_vecs[j];

            std::complex<double> sum{0, 0};
            for (int h = 0; h < nirrep_; h++)
                // Avoid using true_dot until einsums fixes segfault on M-intel runners.
                // sum += einsums::linear_algebra::true_dot(ei[h], ej[h]);
                for (int p = 0; p < irrep_sizes_[h]; p++)
                    for (int q = 0; q < irrep_sizes_[h]; q++)
                        sum += std::conj(ei->block(h).subscript(p, q)) * ej->block(h).subscript(p, q);
            B_(j, i) = sum;
            // B_ is real and symmetric
            if (i != j) B_(i, j) = sum;
        }
    }

    // Container for the column vector solution of gesv. It should be a column vector
    // but has rank == 2. This is a weird artifact of gesv.
    auto C_temp = einsums::Tensor<std::complex<double>, 2>("solution vector", 1, diis_max + 1);

    // Fill the last element of the RHS vector to add the constraint
    C_temp(0, diis_max) = -1;

    // Solves for the DIIS coefficients
    einsums::linear_algebra::gesv(&B_, &C_temp);

    // Extrapolate Fp_ from solution coefficients
    Fp_->zero();
    for (int i = 0; i < diis_max; i++) {
        einsums::linear_algebra::axpy(C_temp(0, i), *Fdiis[i], Fp_.get());
    }
}

void CGHF::form_FDSmSDF() {
    // Computes orbital gradient [F, D] with the overlap matrix EINS_ as the metric by
    // first taking F @ D @ S and storing it in FDSmSDF_
    // S @ D @ F then gets stored in temp2_ and subtracted from FDSmSDF_
    // Then convert to orthogonal basis

    FDSmSDF_->zero();
    temp1_->zero();
    temp2_->zero();

    // F @ D @ S
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, *D_, *EINS_, std::complex<double>{0.0},
                                                temp1_.get());
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, *F_, *temp1_, std::complex<double>{0.0},
                                                FDSmSDF_.get());

    // S @ D @ F
    temp1_->zero();
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, *D_, *F_, std::complex<double>{0.0},
                                                temp1_.get());
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, *EINS_, *temp1_, std::complex<double>{0.0},
                                                temp2_.get());

    (*FDSmSDF_) -= (*temp2_);

    // Orthogonalize gradient [F, D]
    temp1_->zero();
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, *FDSmSDF_, *EINX_, std::complex<double>{0.0},
                                                temp1_.get());
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, *EINX_, *temp1_, std::complex<double>{0.0},
                                               FDSmSDF_.get());

    temp1_->zero();
    temp2_->zero();
}

/* Helper function to return SharedMatrix
 */
SharedMatrix CGHF::get_shared_FDSmSDF() {
    CGHF::form_FDSmSDF();

    Dimension sizes(irrep_sizes_);
    auto out = std::make_shared<Matrix>("FDSmSDF orth", sizes, sizes);

    for (int h = 0; h < sizes.n(); h++) {
        for (int p = 0; p < sizes[h]; p++) {
            for (int q = 0; q < sizes[h]; q++) {
                out->set(h, p, q, (FDSmSDF_->block(h).subscript(p, q)).real());
            }
        }
    }

    return out;
}


/*
/ Equation for GHF S^2 expectation value taken from Nascimento, DePrince
/ https://pubmed.ncbi.nlm.nih.gov/29595979/
/ DOI:  10.1021/acs.jctc.7b01288 
/ 
/ <S^2> = (3/4)N + (1/4)[Tr(D_AA) - Tr(D_BB)]**2 
/       - (1/4)Tr(D_AA @ D_AA) - (1/4)Tr(D_BB @ D_BB)
        + (3/2)Tr(D_AB @ D_BA) - Tr(D_AA @ D_BB)
/
/ where N = Tr(D_AA + D_BB)
*/
void CGHF::s2_expectation_value() {
    double S2 = 0.0;

    for (int h = 0; h < nirrep_; h++) {
        // Gather D_AA, D_BB, D_AB, and D_BA
        einsums::TensorView D_AA = D_->block(h)(einsums::Range{0, nsopi_[h]}, einsums::Range{0, nsopi_[h]});
        einsums::TensorView D_BB = D_->block(h)(einsums::Range{nsopi_[h], irrep_sizes_[h]}, einsums::Range{nsopi_[h], irrep_sizes_[h]});
        einsums::TensorView D_AB = D_->block(h)(einsums::Range{0, nsopi_[h]}, einsums::Range{nsopi_[h], irrep_sizes_[h]});
        einsums::TensorView D_BA = D_->block(h)(einsums::Range{nsopi_[h], irrep_sizes_[h]}, einsums::Range{0, nsopi_[h]});

        double N = 0.0;

        for (int p = 0; p < nsopi_[h]; p++) {
            N += D_AA(p, p).real() + D_BB(p, p).real();
        }

        // Create buffer for matrix-matrix multiplications
        auto mat_buffer = einsums::Tensor<std::complex<double>, 2>("Matrix-matrix container", nsopi_[h], nsopi_[h]);
        mat_buffer.zero();

        //Term 1
        S2 += (3/4) * N;

        //Term 2
        double trAA_minus_trBB = 0.0;
        for (int p = 0; p < nsopi_[h]; p++) {
            trAA_minus_trBB += (D_AA(p, p) - D_BB(p, p)).real();
        }
        S2 += (1/4) * std::pow(trAA_minus_trBB, 2);

        //Term 3: mat_buffer holds D_AA @ D_AA
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, D_AA, D_AA, std::complex<double>{0.0},
                                                &mat_buffer);
        
        double tr_AAAA = 0.0;
        for (int p = 0; p < nsopi_[h]; p++) {
            tr_AAAA += mat_buffer(p, p).real();
        }
        S2 -= (1/4) * tr_AAAA;

        mat_buffer.zero();

        //Term 4: mat_buffer holds D_BB @ D_BB
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, D_BB, D_BB, std::complex<double>{0.0},
                                                &mat_buffer);

        double tr_BBBB = 0.0;
        for (int p = 0; p < nsopi_[h]; p++) {
            tr_BBBB += mat_buffer(p, p).real();
        }
        S2 -= (1/4) * tr_BBBB;

        mat_buffer.zero();

        //Term 5: mat_buffer holds D_AB @ D_BA
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, D_AB, D_BA, std::complex<double>{0.0},
                                                &mat_buffer);

        double tr_ABBA = 0.0;
        for (int p = 0; p < nsopi_[h]; p++) {
            tr_ABBA += mat_buffer(p, p).real();
        }
        S2 += (3/2) * tr_ABBA;

        mat_buffer.zero();

        //Term 6: mat_buffer holds D_AA @ D_BB
        einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, D_AA, D_BB, std::complex<double>{0.0},
                                                &mat_buffer);

        double tr_AABB = 0.0;
        for (int p = 0; p < nsopi_[h]; p++) {
            tr_AABB += mat_buffer(p, p).real();
        }
        S2 -= tr_AABB;
    }

    outfile->Printf("\n\t<S^2> = %.10f\n", S2);

}

void CGHF::finalize() {
    s2_expectation_value();
    HF::finalize();
}

std::shared_ptr<CGHF> CGHF::c1_deep_copy(std::shared_ptr<BasisSet> basis) {
    auto wfn = Wavefunction::c1_deep_copy(basis);

    auto hf_wfn = std::shared_ptr<CGHF>(new CGHF(wfn, functional_, wfn->options(), wfn->psio()));

    return hf_wfn;
}

#endif

}  // namespace scf
}  // namespace psi
