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

#ifdef USING_Einsums

#include "psi4/psi4-dec.h"
#include "psi4/physconst.h"

#include "psi4/libmints/sox2c1e.h"

#include "psi4/libmints/rel_potential.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/libmints/basisset.h"

#include <Einsums/LinearAlgebra.hpp>
#include <Einsums/TensorAlgebra.hpp>
#include <Einsums/Runtime.hpp>
#include <Einsums/Concepts/TensorConcepts.hpp>

#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

namespace {

/* Computes the conjugate transpose of const reference A and stores result in B */
void conj_T(ComplexMatrix const& A, std::shared_ptr<ComplexMatrix> B) {
    // It may seem kitsch to pass a dereferenced shared_ptr to `A`, but it's safe.
    // `A` can neither be changed, nor become a dangling reference.
    for (int h = 0; h < A.num_blocks(); h++) {
        einsums::tensor_algebra::permute<true>(
            std::complex<double>{0.0}, einsums::Indices{einsums::index::i, einsums::index::j}, &B->block(h),
            std::complex<double>{1.0}, einsums::Indices{einsums::index::j, einsums::index::i}, A.block(h));
    }
}

/* Convenient helper functions for gemm call of complex underlying type.
 * Templating allows for the expected
 *   gemm(ComplexMatrix const& A, ComplexMatrix const& B,
 *        std::shared_ptr<ComplexMatrix> C);
 *
 * as well as intermingling TensorView types:
 *   gemm(TensorView<ComplexMatrix> const& A, TensorView<ComplexMatrix> const& B,
 *        ComplexMatrix* C);
 *
 * and shared pointers
 *   gemm(TensorView<ComplexMatrix> const& A, TensorView<ComplexMatrix> const& B,
 *        std::shared_ptr<ComplexMatrix> C);
 */
template <einsums::TensorConcept T1, einsums::TensorConcept T2>
    requires einsums::SameUnderlying<einsums::RemoveViewT<T1>, T2>
void gemm(T1 const& A, T1 const& B, T2* C) {
    einsums::linear_algebra::gemm<false, false>(std::complex<double>(1.0), A, B, std::complex<double>(0.0), C);
}

template <einsums::TensorConcept T1, typename T2>
    requires einsums::SameUnderlying<einsums::RemoveViewT<T1>, T2>
void gemm(T1 const& A, T1 const& B, std::shared_ptr<T2> C) {
    einsums::linear_algebra::gemm<false, false>(std::complex<double>(1.0), A, B, std::complex<double>(0.0), C.get());
}

}  // namespace

namespace psi {

SOX2C1e::SOX2C1e(std::shared_ptr<BasisSet> basis, std::shared_ptr<BasisSet> x2c_basis)
    : contractedBasis_(basis),  // BASIS
      uncBasis_(x2c_basis)      // BASIS_RELATIVISTIC
{
    outfile->Printf(R"(
      _____ ____ _  _____   ______    ___
     / ___// __ \ |/ /__ \ / ____/   <  /__
     \__ \/ / / /   /__/ // /  ______/ / _ \
    ___/ / /_/ /   |/ __// /__/_____/ /  __/
   /____/\____/_/|_/____/\____/    /_/\___/
            by Nathan Gillispie
)");

    outfile->Printf("\n  ==> X2C Options <==");
    outfile->Printf("\n    Ref Basis: " + contractedBasis_->name());
    outfile->Printf("\n    X2C Basis: " + uncBasis_->name());
    outfile->Printf("\n    The X2C Hamiltonian will be computed in the X2C Basis.\n\n");

    // Turn off non-critical einsums logs to stdout
    std::vector<std::string> ein_argv{"psi4", "--einsums:no-profiler-report", "--einsums:log-level", "3"};
    einsums::initialize(ein_argv);
}

SOX2C1e::~SOX2C1e() {}

void SOX2C1e::compute(SharedMatrix S, SharedMatrix T, SharedMatrix V, SharedMatrix Hx, SharedMatrix Hy,
                      SharedMatrix Hz) {
    outfile->Printf("\033[41m   1. Compute integrals\033[0m\n");
    // TODO: find out why we are using the uncontracted basis to get the SO integrals??
    auto integral = std::make_shared<IntegralFactory>(uncBasis_, uncBasis_, uncBasis_, uncBasis_);

    // Obtain the dimension object to initialize the factory.
    auto soBasis = std::make_shared<SOBasisSet>(uncBasis_, integral);
    nsopi_ = soBasis->dimension();
    nsopi2c_ = nsopi_ + nsopi_;
    nsopi4c_ = nsopi2c_ + nsopi2c_;

    auto soFactory = std::make_shared<MatrixFactory>();
    soFactory->init_with(soBasis);
    compute_integrals(integral, soFactory);

    outfile->Printf("\033[41m   2. Form Dirac Hamiltonian\033[0m\n");
    auto Hdirac = std::make_shared<ComplexMatrix>("Dirac Hamiltonian", nsopi4c_.blocks());
    form_dirac_hamiltonian(Hdirac);

    outfile->Printf("\033[41m   3. Diagonalize Dirac Hamiltonian\033[0m\n");
    auto orth = std::make_shared<ComplexMatrix>("4c metric^(-1/2)", nsopi4c_.blocks());
    auto orth_H = std::make_shared<ComplexMatrix>("4c metric^(-1/2).conj().T", nsopi4c_.blocks());
    form_orth(orth);
    conj_T(*orth, orth_H);

    auto Hevec = std::make_shared<ComplexMatrix>("Dirac eigenvectors", nsopi4c_.blocks());
    auto Heval = einsums::BlockTensor<double, 1>("Dirac eigenvalues", nsopi4c_.blocks());
    auto tmp = std::make_shared<ComplexMatrix>("tmp", nsopi4c_.blocks());

    // Orthogonalize Hdirac
    gemm(*orth_H, *Hdirac, tmp);
    gemm(*tmp, *orth, Hevec);

    for (int h = 0; h < nsopi_.n(); h++) {
        if (nsopi_[h] == 0) continue;

        einsums::linear_algebra::heev<true>(&Hevec->block(h), &Heval[h]);

        for (int i = 0; i < nsopi4c_[h]; i++) {
            if ((i % 5) == 0) {
                outfile->Printf("\n  (%03d) ", i);
            }
            outfile->Printf("  %13.6f", Heval[h](i));
        }
        outfile->Printf("\n");
    }
    conj_T(*Hevec, tmp);       // necessary due to different conventions in BLAS & Einsums
    gemm(*orth, *tmp, Hevec);  // Back-transform Hevec

    outfile->Printf("\033[41m   4. Form X & S\033[0m\n");
    auto X = std::make_shared<ComplexMatrix>("X", nsopi2c_.blocks());
    form_X(*Hevec, X);

    auto X_H = std::make_shared<ComplexMatrix>("X.conj().T", nsopi2c_.blocks());
    auto X_HT = std::make_shared<ComplexMatrix>("X.conj().T @ T", nsopi2c_.blocks());
    auto TX = std::make_shared<ComplexMatrix>("T @ X", nsopi2c_.blocks());
    auto X_HTX = std::make_shared<ComplexMatrix>("X.conj().T @ T @ X", nsopi2c_.blocks());
    conj_T(*X, X_H);
    gemm(*X_H, *kinetic_, X_HT);
    conj_T(*X_HT, TX);
    gemm(*X_HT, *X, X_HTX);

    auto Stilde = std::make_shared<ComplexMatrix>("Stilde", nsopi2c_.blocks());
    form_Stilde(*X_HTX, Stilde);


    outfile->Printf("\033[41m   5. Form R\033[0m\n");
    auto R = std::make_shared<ComplexMatrix>("R", nsopi2c_.blocks());
    form_R(*Stilde, *orth, R);


    outfile->Printf("\033[41m   6. Form hFW+\033[0m\n");
    auto hFW = std::make_shared<ComplexMatrix>("hFW+", nsopi2c_.blocks());
    (*hFW) = (*TX);
    (*hFW) += (*X_HT);
    (*hFW) -= (*X_HTX);
    (*hFW) += (*nuclear_);

    // Now for 1/4c²X†WX term
    tmp.reset();
    tmp = std::make_shared<ComplexMatrix>("tmp", nsopi2c_.blocks()); // Different size tmp
    gemm(*X_H, *rel_pot_, tmp);
    gemm(*tmp, *X, rel_pot_);
    (*hFW) += (*rel_pot_);

    Hevec.reset();
    Hevec = std::make_shared<ComplexMatrix>(*hFW);

    conj_T(*R, tmp);
    gemm(*tmp, *hFW, Hevec);
    gemm(*Hevec, *R, hFW);

    outfile->Printf("\033[41m   7. Project\033[0m\n");
    outfile->Printf("\033[41m   8. Test hFW+\033[0m\n");
    Heval = einsums::BlockTensor<double, 1>("X2C eigenvalues", nsopi2c_.blocks());

    orth.reset();
    orth = std::make_shared<ComplexMatrix>("none", nsopi2c_.blocks());
    orth->zero();

    for (int h = 0; h < nsopi_.n(); ++h) {
        int dim = nsopi_[h];
        auto& orth_block = orth->block(h);
        for (int p = 0; p < dim; p++) {
            for (int q = 0; q < dim; q++) {
                orth_block(p, q) = sMat->get(h, p, q);
                orth_block(p + dim, q + dim) = sMat->get(h, p, q);
            }
        }
    }

    orth_H.reset();
    orth_H = std::make_shared<ComplexMatrix>("orthH", nsopi2c_.blocks());
    conj_T(*orth, orth_H);
    gemm(*orth_H, *hFW, tmp);
    gemm(*tmp, *orth, Hevec);

    for (int h = 0; h < nsopi_.n(); h++) {
        if (nsopi_[h] == 0) continue;

        einsums::linear_algebra::heev<true>(&Hevec->block(h), &Heval[h]);

        for (int i = 0; i < nsopi2c_[h]; i++) {
            if ((i % 5) == 0) {
                outfile->Printf("\n  (%03d) ", i);
            }
            outfile->Printf("  %13.6f", Heval[h](i));
        }
        outfile->Printf("\n");
    }

    outfile->Printf("\033[41m   9. Make real, set to original matrices\033[0m\n");
}

/*
 * Use Libint2 to compute integrals. Store pauli-decomposed relativistic potential W
 * as 4 real SharedMatrix objects.
 */
void SOX2C1e::compute_integrals(std::shared_ptr<IntegralFactory> integral, std::shared_ptr<MatrixFactory> soFactory) {
    std::unique_ptr<OneBodySOInt> sOBI(integral->so_overlap());
    std::unique_ptr<OneBodySOInt> tOBI(integral->so_kinetic());
    std::unique_ptr<OneBodySOInt> vOBI(integral->so_potential());
    std::unique_ptr<OneBodySOInt> wOBI(integral->so_rel_potential());

    // Keeping sMat, tMat for later (computing metric/orth)
    sMat = soFactory->create_shared_matrix("Overlap");
    tMat = soFactory->create_shared_matrix("Kinetic");
    auto vMat = soFactory->create_shared_matrix("Potential");
    auto w0Mat = soFactory->create_shared_matrix("Relativistic Potential (scalar)");
    auto wxMat = soFactory->create_shared_matrix("Relativistic Potential (pauli-x)");
    auto wyMat = soFactory->create_shared_matrix("Relativistic Potential (pauli-y)");
    auto wzMat = soFactory->create_shared_matrix("Relativistic Potential (pauli-z)");

    std::vector<SharedMatrix> wMats{w0Mat, wxMat, wyMat, wzMat};

    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);
    wOBI->compute(wMats);

    const double C2 = .25 / (pc_c_au * pc_c_au);
    for (auto& w : wMats) {
        w->scale(C2);
    }

    // Matrix to einsums conversion: block diag for most except rel_pot_.
    overlap_ = std::make_shared<ComplexMatrix>("Overlap 2c", nsopi2c_.blocks());
    kinetic_ = std::make_shared<ComplexMatrix>("Kinetic 2c", nsopi2c_.blocks());
    nuclear_ = std::make_shared<ComplexMatrix>("Nuclear 2c", nsopi2c_.blocks());
    rel_pot_ = std::make_shared<ComplexMatrix>("s.pVs.p 2c", nsopi2c_.blocks());

    overlap_->zero();
    kinetic_->zero();
    nuclear_->zero();
    rel_pot_->zero();

    const std::complex<double> i{0, 1};
    for (int h = 0; h < nsopi2c_.n(); h++) {
        int nc = nsopi_[h];  // size of 1-component
        auto& Sb = overlap_->block(h);
        auto& Tb = kinetic_->block(h);
        auto& Vb = nuclear_->block(h);
        auto& Wb = rel_pot_->block(h);
        for (int p = 0; p < nc; p++) {
            for (int q = 0; q < nc; q++) {
                Sb(p, q) += sMat->get(h, p, q);
                Sb(p + nc, q + nc) += sMat->get(h, p, q);

                Tb(p, q) += tMat->get(h, p, q);
                Tb(p + nc, q + nc) += tMat->get(h, p, q);

                Vb(p, q) += vMat->get(h, p, q);
                Vb(p + nc, q + nc) += vMat->get(h, p, q);

                Wb(p, q) += wMats[0]->get(h, p, q) + i * wMats[3]->get(h, p, q);
                Wb(p + nc, q + nc) += wMats[0]->get(h, p, q) - i * wMats[3]->get(h, p, q);
                Wb(p, q + nc) += wMats[2]->get(h, p, q) + i * wMats[1]->get(h, p, q);
                Wb(p + nc, q) -= wMats[2]->get(h, p, q) - i * wMats[1]->get(h, p, q);
            }
        }
    }
}

/*   Form the Dirac Hamiltonian
 *        | V       T       |
 *    d = |                 |
 *        | T  1/4c^2 W - T |
 */
void SOX2C1e::form_dirac_hamiltonian(std::shared_ptr<ComplexMatrix> Hdirac) {
    Hdirac->zero();

    auto tmp = std::make_shared<ComplexMatrix>(*rel_pot_);
    // recall that rel_pot_ has already been scaled by .25 / (pc_c_au * pc_c_au)
    (*tmp) -= (*kinetic_);

    for (int h = 0; h < nsopi_.n(); h++) {
        int nc = nsopi2c_[h];
        auto& Hblock = Hdirac->block(h);

        einsums::Range a{0, nsopi2c_[h]};
        einsums::Range b{nsopi2c_[h], nsopi4c_[h]};

        Hblock(a, a) = nuclear_->block(h);
        Hblock(a, b) = kinetic_->block(h);
        Hblock(b, a) = kinetic_->block(h);
        Hblock(b, b) = tmp->block(h);
    }
}

/* Form metric^{-1/2}:
 *    | S        0 |^-1/2
 *    |            |
 *    | 0  T/2c**2 |
 */
void SOX2C1e::form_orth(std::shared_ptr<ComplexMatrix> orth) {
    orth->zero();

    tMat->scale(.5 / (pc_c_au * pc_c_au));
    sMat->power(-.5);
    tMat->power(-.5);

    for (int h = 0; h < nsopi_.n(); ++h) {
        int dim = nsopi_[h];
        auto& orth_block = orth->block(h);
        for (int p = 0; p < dim; p++) {
            for (int q = 0; q < dim; q++) {
                orth_block(p, q) = sMat->get(h, p, q);
                orth_block(p + dim, q + dim) = sMat->get(h, p, q);
                orth_block(p + 2 * dim, q + 2 * dim) = tMat->get(h, p, q);
                orth_block(p + 3 * dim, q + 3 * dim) = tMat->get(h, p, q);
            }
        }
    }
}

void SOX2C1e::form_X(ComplexMatrix const& Hevec, std::shared_ptr<ComplexMatrix> X) {
    for (int h = 0; h < nsopi_.n(); h++) {
        einsums::Range small{2 * nsopi_[h], 4 * nsopi_[h]};
        auto pe = small;  // Positive energy solutions
        einsums::Range large{0, 2 * nsopi_[h]};

        einsums::TensorView C_large = Hevec[h](large, pe);
        einsums::TensorView C_small = Hevec[h](small, pe);

        einsums::linear_algebra::invert(&C_large);
        gemm(C_small, C_large, &X->block(h));
    }
}

/* Stilde = S + \frac{1}{2c^2}X^\dagger T X */
void SOX2C1e::form_Stilde(ComplexMatrix const& X_HTX, std::shared_ptr<ComplexMatrix> Stilde) {
    const std::complex<double> C2{.5 / (pc_c_au * pc_c_au)};

    Stilde->zero();
    Stilde->operator+=(X_HTX);
    einsums::linear_algebra::scale(C2, Stilde.get());

    Stilde->operator+=(*overlap_);
}

/* Form R = S^{-1/2}[S^{-1/2} \tilde S S^{-1/2}]^{-1/2} S^{1/2} */
void SOX2C1e::form_R(ComplexMatrix const& Stilde, ComplexMatrix const& orth, std::shared_ptr<ComplexMatrix> R) {
    // Shalf == S^{1/2}
    auto Shalf = std::make_shared<ComplexMatrix>(orth);
    einsums::linear_algebra::invert(Shalf.get());

    auto tmp = std::make_shared<ComplexMatrix>(Stilde);
    auto evec = std::make_shared<ComplexMatrix>(Stilde);
    auto eval = einsums::BlockTensor<double, 1>("eigenvalues", nsopi2c_.blocks());

    // tmp <- [S^{-1/2} \tilde S S^{-1/2}]
    gemm(orth, Stilde, evec);
    gemm(*evec, orth, tmp);

    // Computing [S^{-1/2} \tilde S S^{-1/2}]^{-1/2}

    // TODO: when einsums pow bug is fixed, switch to this:
    // (*R) = einsums::linear_algebra::pow(*tmp, -.5);
    // For now, just do it the old-fashioned way.

    for (int h = 0; h < nsopi_.n(); h++) {
        if (nsopi_[h] == 0) continue;
        einsums::linear_algebra::heev<true>(&tmp->block(h), &eval[h]);
    }
    conj_T(*tmp, evec);

    // R <- Σ Λ^{-1/2} Σ^\dagger == evec @ np.diag(eval) @ tmp
    // scale each column in evec by eval: evec @ np.diag(eval)
    for (int h = 0; h < nsopi2c_.n(); h++) {
        auto& vecblock = evec->block(h);
        for (int i = 0; i < nsopi2c_[h]; i++) {
            einsums::linear_algebra::scale_column(i, std::pow(eval[h](i), -.5), &vecblock);
        }
    }
    gemm(*evec, *tmp, R);

    // S^{-1/2}...S^{1/2}
    gemm(orth, *R, tmp);
    gemm(*tmp, *Shalf, R);
}

}  // namespace psi

#endif  // USING_Einsums
