/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "cghf.h"

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/orthog.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libfock/ComplexJK.h"
#include "psi4/libfock/v.h"

#include <Einsums/Print.hpp>
// DEBUG STATEMENTS
#include <pybind11/pybind11.h>
// END DEBUG

namespace {

// Takes an nsopi_-shaped square SharedMatrix and copies to 2 (two) diagonal
// blocks **per irrep** into each tile of the provided ComplexMatrix.
void copy_matrix_to_complex(const psi::Matrix& A, psi::ComplexMatrix& B) {
    const int nirrep = A.nirrep();
    std::vector<int> row_dim(nirrep);
    std::vector<int> col_dim(nirrep);

    for (int h = 0; h < nirrep; h++) {
        row_dim[h] = A.rowspi(h) * 2;
        col_dim[h] = A.colspi(h) * 2;
    }

    B = psi::ComplexMatrix{B.name(), row_dim, col_dim};
    B.zero();

    for (int h = 0; h < nirrep; h++) {
        const int r_ = A.rowspi(h);
        const int c_ = A.colspi(h);
        for (int i = 0; i < r_; i++) {
            for (int j = 0; j < c_; j++) {
                B.tile(h, h)(i, j) = A(h, i, j);
                B.tile(h, h)(i + r_, j + c_) = A(h, i, j);
            }
        }
    }
}

}

namespace psi {
namespace scf {

CGHF::CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> func)
    : CGHF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {}

CGHF::CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
           std::shared_ptr<PSIO> psio)
    : ComplexWavefunction(options), BaseHF(func) {
    // Analogous to HF::HF's shallow_copy(ref_wfn); full ComplexWavefunction::shallow_copy comes later.
    molecule_ = ref_wfn->molecule();
    basisset_ = ref_wfn->basisset();
    psio_ = psio;
    // Options-only base ctor skipped init; run it now that molecule/basis are set.
    ComplexWavefunction::common_init();
    common_init();
}

CGHF::~CGHF() {}

void CGHF::common_init() {
    BaseHF::common_init(options_, module_, molecule_, dipole_field_strength_);
    name_ = "CGHF";
    // Prefer BaseHF::nelectron_ over BaseWavefunction's; copy from ComplexWavefunction.
    nelectron_ = nelec_;

    if (molecule_->schoenflies_symbol() != "c1" || nirrep_ != 1) {
        throw PSIEXCEPTION("CGHF currently supports only C1 symmetry. Set symmetry c1 in the molecule block.");
    }

    // DFT stuff (would typically go in subclass_init)
    setup_potential();

    std::vector<size_t> irrep_sizes(nirrep_);
    for (int h = 0; h < nirrep_; h++) {
        irrep_sizes[h] = static_cast<size_t>(nsopi_[h] * 2);
        // nelecpi_[h] = static_cast<size_t>(nalphapi_[h] + nbetapi_[h]);
        nelecpi_[0] = static_cast<size_t>(nelectron_);
        if (h > 0) throw PSIEXCEPTION("looking for something? *waves nalphapi_ in front of face*");
    }

    T_ = std::make_shared<ComplexMatrix>("T", irrep_sizes, irrep_sizes);
    V_ = std::make_shared<ComplexMatrix>("V", irrep_sizes, irrep_sizes);
    H_ = std::make_shared<ComplexMatrix>("H", irrep_sizes, irrep_sizes);
    F_ = std::make_shared<ComplexMatrix>("F", irrep_sizes, irrep_sizes);
    G_ = std::make_shared<ComplexMatrix>("G", irrep_sizes, irrep_sizes);
    J_ = std::make_shared<ComplexMatrix>("J", irrep_sizes, irrep_sizes);
    K_ = std::make_shared<ComplexMatrix>("K", irrep_sizes, irrep_sizes);

    T_->zero();
    V_->zero();
    H_->zero();
    F_->zero();
    G_->zero();
    J_->zero();
    K_->zero();

    // We don't know the sizes of these until nmopi_ fills in form_Shalf();
    S_ = std::make_shared<ComplexMatrix>(); S_->set_name("Overlap");
    X_ = std::make_shared<ComplexMatrix>(); X_->set_name("Orthogonalization");
    // C_ is resized in form_C once X_ is known; D_ matches SO dimension.
    C_ = std::make_shared<ComplexMatrix>(); C_->set_name("MO coefficients");
    D_ = std::make_shared<ComplexMatrix>("SCF density", irrep_sizes, irrep_sizes);
    D_->zero();
}

void CGHF::setup_potential() {
    if (functional_->needs_xc()) {
        throw PSIEXCEPTION("DFT functionals are not supported with CGHF.");
    } else {
        potential_ = nullptr;
    }
}

void CGHF::set_jk(std::shared_ptr<ComplexJK> jk) {
    // Cheap basis check
    int jk_nbf = jk->basisset()->nbf();
    int hf_nbf = basisset_->nbf();
    if (hf_nbf != jk_nbf) {
        throw PSIEXCEPTION("Tried setting a JK object whos number of basis functions does not match HF's!");
    }

    jk_ = jk;
}

void CGHF::form_H() {
    SharedMatrix T_real = mintshelper()->so_kinetic();
    SharedMatrix V_real = mintshelper()->so_potential();

    copy_matrix_to_complex(*T_real, *T_);
    copy_matrix_to_complex(*V_real, *V_);

    if (debug_ > 2) einsums::fprintln(*outfile->stream(), *T_);
    if (debug_ > 2) einsums::fprintln(*outfile->stream(), *V_);

    if (perturb_h_)
        throw PSIEXCEPTION("Perturbed Hamiltonian not supported with CGHF.");

    if (external_pot_)
        throw PSIEXCEPTION("External potential not supported with CGHF.");

    (*H_) = (*T_);
    (*H_) += (*V_);

    if (print_ > 3) einsums::fprintln(*outfile->stream(), *H_);
}

void CGHF::form_Shalf() {
    BasisSetOrthogonalization::OrthogonalizationMethod method;
    if (options_.get_str("S_ORTHOGONALIZATION") == "SYMMETRIC")
        method = BasisSetOrthogonalization::Symmetric;
    else if (options_.get_str("S_ORTHOGONALIZATION") == "CANONICAL")
        method = BasisSetOrthogonalization::Canonical;
    else if (options_.get_str("S_ORTHOGONALIZATION") == "PARTIALCHOLESKY")
        method = BasisSetOrthogonalization::PartialCholesky;
    else if (options_.get_str("S_ORTHOGONALIZATION") == "AUTO")
        method = BasisSetOrthogonalization::Automatic;
    else
        throw PSIEXCEPTION("Unrecognized S_ORTHOGONALIZATION method\n");


    double lindep_tolerance = options_.get_double("S_TOLERANCE");
    double cholesky_tolerance = options_.get_double("S_CHOLESKY_TOLERANCE");

    SharedMatrix S_temp = mintshelper()->so_overlap();
    BasisSetOrthogonalization orthog(method, S_temp, lindep_tolerance, cholesky_tolerance, print_);

    // Transform
    SharedMatrix X_temp = orthog.basis_to_orthog_basis();

    // Update nmo_
    auto nmopi_ = X_temp->colspi();
    auto nmo_ = nmopi_.sum();

    // Double check occupation vectors
    for (int h = 0; h < X_temp->nirrep(); ++h) {
        if (nelecpi_[h] > nmopi_[h]) {
            throw PSIEXCEPTION("Not enough molecular orbitals to satisfy requested occupancies");
        }
    }

    // Create the correct sized quantities now.
    copy_matrix_to_complex(*S_temp, *S_);
    copy_matrix_to_complex(*X_temp, *X_);
}

void CGHF::guess() {
    double guess_E;
    std::string guess_type = options_.get_str("GUESS");

    if (guess_type == "CORE") {
        if (print_) outfile->Printf("  SCF Guess: Core (One-Electron) Hamiltonian.\n\n");

        (*F_) = (*H_); // Try the core Hamiltonian as the Fock Matrix
        form_initial_C(); // calls only CGHF::form_C()
        form_D();
        guess_E = compute_initial_E();
    } else {
        throw PSIEXCEPTION("CGHF '" + guess_type + "' GUESS not implemented. Use 'CORE'.");
    }

    energies_["Total Energy"] = 0.0;  // don't use this guess in our convergence checks
}

double CGHF::compute_initial_E() {
    std::complex<double> E;

    /*  \sum_{ij}h_{ij} \gamma_{ji} */
    using namespace einsums;
    tensor_algebra::einsum(Indices{}, &E, Indices{index::i, index::j}, *H_, Indices{index::j, index::i}, *D_);

    if (E.imag() > 1e-12) {
        outfile->Printf("WARNING: CGHF::compute_initial_E found large imaginary %+fi\n"
                        "  Is D Hermitian?\n", E.imag());
    }

    return nuclearrep_ + E.real();
}

void CGHF::form_C(double shift) {
    if (shift != 0.0) throw PSIEXCEPTION("Level shifting not available for CGHF.");

    auto temp = ComplexMatrix("temp", F_->tile_size(0), X_->tile_size(1));
    auto XFX = ComplexMatrix("Othogonalized Fock", X_->tile_size(1), X_->tile_size(1));

    // using namespace einsums;

    // Form F' = X'FX for canonical orthogonalization
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, *F_, *X_,
                                       std::complex<double>{0.0}, &temp);
    einsums::linear_algebra::gemm<true, false>(std::complex<double>{1.0}, *X_, temp,
                                      std::complex<double>{0.0}, &XFX);

    // Form C' = eig(F')
    temp = ComplexMatrix("temp", XFX.tile_sizes());
    temp.zero();
    for (int h = 0; h < nirrep_; h++) {
        // Do not diagonalize 0x0 matrix
        if (nsopi_[h] == 0) continue;

        auto evals = einsums::Tensor<double, 1>("Fock evals", nsopi_[h]*2);
        evals.zero();

        // Hermitian eigensolver one block at a time
        einsums::linear_algebra::heev<true>(&XFX.tile(h, h), &evals);

        // TODO: initialize epsilon_ with correct shape then copy evals->*epsilon_

        // heev retuns the wrong side, so we need to take the inverse (hermitian adjoint)

        // Takes the conjugate transpose of XFX (e.g. ij -> ji) to give us the proper eigenvectors
        // NOTE: the template parameters <true> states to take the conjugate (Einsums v1.x)
        einsums::tensor_algebra::permute<true>(
            std::complex<double>{0.0}, einsums::Indices{einsums::index::i, einsums::index::j}, &temp.tile(h,h),
            std::complex<double>{1.0}, einsums::Indices{einsums::index::j, einsums::index::i}, XFX.tile(h,h));
    }

    // Form C_ := X_ @ temp
    C_ = std::make_shared<ComplexMatrix>("MO coefficients", X_->tile_sizes());
    einsums::linear_algebra::gemm<false, false>(std::complex<double>{1.0}, *X_, temp,
                                       std::complex<double>{0.0}, C_.get());

    // find_occupation();
}

void CGHF::form_D() {
    throw pybind11::attribute_error("form_D is not implemented for CGHF.");
}

void CGHF::form_F() {
    throw pybind11::attribute_error("form_F is not implemented for CGHF.");
}

void CGHF::form_G() {
    throw pybind11::attribute_error("form_G is not implemented for CGHF.");
}

std::shared_ptr<ComplexJK> CGHF::build_jk(size_t memory) const {
    auto jk = ComplexJK::build_JK(get_basisset("ORBITAL"), get_basisset("DF_BASIS_SCF"),
                                  Process::environment.options, false, memory);
    return jk;
}

}  // namespace scf
}  // namespace psi
