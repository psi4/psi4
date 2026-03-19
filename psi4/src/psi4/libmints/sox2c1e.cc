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

#include "psi4/libpsi4util/process.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

/* Computes the conjugate transpose of const reference A and stores result in B */
void conj_T(ComplexMatrix const& A, SharedComplexMatrix B) {
    // It may seem kitsch to pass a dereferenced shared_ptr to `A`, but it's safe.
    // `A` can neither be changed, nor become a dangling reference.
    // TODO: permute BlockTensor directly when Einsums bug is fixed
    for (int h = 0; h < A.num_blocks(); h++) {
        if (A.block_dims()[h] == 0) continue;
        einsums::tensor_algebra::permute<true>(
            std::complex<double>{0.0}, einsums::Indices{einsums::index::i, einsums::index::j}, &B->block(h),
            std::complex<double>{1.0}, einsums::Indices{einsums::index::j, einsums::index::i}, A.block(h));
    }
}

/* Convenient helper functions for gemm call of complex underlying type.
 * Templating allows for the expected
 *   gemm(ComplexMatrix const& A, ComplexMatrix const& B,
 *        SharedComplexMatrix C);
 *
 * as well as intermingling TensorView types:
 *   gemm(TensorView<ComplexMatrix> const& A, TensorView<ComplexMatrix> const& B,
 *        ComplexMatrix* C);
 *
 * and shared pointers
 *   gemm(TensorView<ComplexMatrix> const& A, TensorView<ComplexMatrix> const& B,
 *        SharedComplexMatrix C);
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
    : aoBasis_(basis),        // BASIS
      deconBasis_(x2c_basis)  // BASIS_RELATIVISTIC
{
    outfile->Printf(R"(
     _____ ____ _  _____   ______    ___
    / ___// __ \ |/ /__ \ / ____/   <  /__
    \__ \/ / / /   /__/ // /  ______/ / _ \
   ___/ / /_/ /   |/ __// /__/_____/ /  __/
  /____/\____/_/|_/____/\____/    /_/\___/
           by Nathan Gillispie
)");

    outfile->Printf("\n   ==> SOX2C-1e Options <==");
    outfile->Printf("\n     Ref Basis: " + aoBasis_->name());
    outfile->Printf("\n     X2C Basis: " + deconBasis_->name());
    outfile->Printf("\n     The X2C Hamiltonian will be computed in the X2C Basis.\n\n");

    // Turn off non-critical einsums logs to stdout
    std::vector<std::string> ein_argv{"psi4", "--einsums:no-profiler-report", "--einsums:log-level", "3"};
    einsums::initialize(ein_argv);
}

SOX2C1e::~SOX2C1e() {}

/*
 * Main function called in mintshelper.cc. S, T, V, \vec{H} are set here.
 */
void SOX2C1e::compute(SharedMatrix S, SharedMatrix T, SharedMatrix V, SharedMatrix Hx, SharedMatrix Hy,
                      SharedMatrix Hz) {
    /* 1. Compute integrals. */
    auto integral = std::make_shared<IntegralFactory>(deconBasis_, deconBasis_, deconBasis_, deconBasis_);

    auto deconSoBasis = std::make_shared<SOBasisSet>(deconBasis_, integral);
    nsopi_ = deconSoBasis->dimension();
    nsopi2c_ = nsopi_ + nsopi_;
    nsopi4c_ = nsopi2c_ + nsopi2c_;

    auto soFactory = std::make_shared<MatrixFactory>();
    soFactory->init_with(deconSoBasis);
    compute_integrals(integral, soFactory);

    /* 2. Form 1e Dirac Hamiltonian. */
    auto Hdirac = std::make_shared<ComplexMatrix>("Dirac Hamiltonian", nsopi4c_.blocks());
    form_dirac_hamiltonian(Hdirac);

    /* 3. Diagonalize 1e Dirac Hamiltonian with metric */
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

        if (Process::environment.options.get_int("PRINT") > 1) {
            outfile->Printf("  Eigenvalues of Dirac Hamiltonian. Irrep %d", h+1);
            for (int i = 0; i < nsopi4c_[h]; i++) {
                if ((i % 5) == 0) {
                    outfile->Printf("\n  (%03d) ", i);
                }
                outfile->Printf("  %13.6f", Heval[h](i));
            }
            outfile->Printf("\n\n");
        }
    }
    // Must do conj_T due to different conventions in BLAS & Einsums
    conj_T(*Hevec, tmp);
    gemm(*orth, *tmp, Hevec);  // Back-transform
    tmp.reset(); // Last time using 4c matrices

    /* 4. Form X and Stilde */
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

    /* 5. Form R */
    auto R = std::make_shared<ComplexMatrix>("R", nsopi2c_.blocks());
    form_R(*Stilde, *orth, R);

    /* 6. Form NESC FW+ Hamiltonian (Split into V and T) */
    auto Vx2c = ComplexMatrix("Vx2c", nsopi2c_.blocks());
    auto Tx2c = ComplexMatrix("Tx2c", nsopi2c_.blocks());
    Tx2c = (*TX);
    Tx2c += (*X_HT);
    Tx2c -= (*X_HTX);

    Vx2c = (*nuclear_);

    // Now for 1/4c²X†WX term
    tmp = std::make_shared<ComplexMatrix>("tmp", nsopi2c_.blocks());  // Different size tmp
    tmp->zero();
    gemm(*X_H, *rel_pot_, tmp);
    gemm(*tmp, *X, rel_pot_);
    Vx2c += (*rel_pot_);

    // F = R† hFW R
    // X_H is now R.conj().T
    conj_T(*R, X_H);

    gemm(*X_H, Vx2c, tmp);
    gemm(*tmp, *R, &Vx2c);

    gemm(*X_H, Tx2c, tmp);
    gemm(*tmp, *R, &Tx2c);

    /* 7. Express Tx2c + Vx2c as real matrices */
    // Decontracted basis SharedMatrices
    auto S_x2c = std::make_shared<Matrix>(sMat->clone());
    auto T_x2c = std::make_shared<Matrix>("Tx2c", nsopi_, nsopi_);
    auto V_x2c = std::make_shared<Matrix>("Vx2c", nsopi_, nsopi_);
    auto Hx_x2c = std::make_shared<Matrix>("Hx", nsopi_, nsopi_);
    auto Hy_x2c = std::make_shared<Matrix>("Hy", nsopi_, nsopi_);
    auto Hz_x2c = std::make_shared<Matrix>("Hz", nsopi_, nsopi_);
    form_pauli(Tx2c, Vx2c, T_x2c, V_x2c, Hx_x2c, Hy_x2c, Hz_x2c);

    /* 8. Project? */
    if (aoBasis_->name() != deconBasis_->name()) {
        SharedMatrix D = get_projection();
        S_x2c->transform(D);
        T_x2c->transform(D);
        V_x2c->transform(D);
        Hx_x2c->transform(D);
        Hy_x2c->transform(D);
        Hz_x2c->transform(D);
    } else {
        outfile->Printf("\n  NOTE: Not projecting RELATIVISTIC_BASIS X2C integrals"
                        " to BASIS. Basis sets appear to be the same.\n");
    }

    S->copy(S_x2c);
    T->copy(T_x2c);
    V->copy(V_x2c);
    Hx->copy(Hx_x2c);
    Hy->copy(Hy_x2c);
    Hz->copy(Hz_x2c);

    /* 9. Test hFW+ */
    bool same_evals = test_hFW(Heval, S, T, V, Hx, Hy, Hz);

    if (!same_evals) {
        outfile->Printf("\n  WARNING: The X2C and Dirac Hamiltonians have substatially different eigenvalues!\n");
        outfile->Printf("           This is probably caused by the recontraction of the basis set.\n\n");
    }

    outfile->Printf("  ==> SOX2C-1e Integrals computed successfully! <==\n\n");
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

    // W always appears with a prefactor when it is used here.
    const double C2 = .25 / (pc_c_au * pc_c_au);
    for (auto& w : wMats) {
        w->scale(C2);
    }

    // Matrix to einsums conversion: block diag for most except rel_pot_.
    overlap_ = std::make_shared<ComplexMatrix>("Overlap 2c", nsopi2c_.blocks());
    kinetic_ = std::make_shared<ComplexMatrix>("Kinetic 2c", nsopi2c_.blocks());
    nuclear_ = std::make_shared<ComplexMatrix>("Nuclear 2c", nsopi2c_.blocks());
    rel_pot_ = std::make_shared<ComplexMatrix>("s.pVs.p 2c", nsopi2c_.blocks());

    const std::complex<double> i{0, 1};
    for (int h = 0; h < nsopi_.n(); h++) {
        int nc = nsopi_[h];  // size of 1-component
        auto& Sb = overlap_->block(h);
        auto& Tb = kinetic_->block(h);
        auto& Vb = nuclear_->block(h);
        auto& Wb = rel_pot_->block(h);
        for (int p = 0; p < nc; p++) {
            for (int q = 0; q < nc; q++) {
                Sb(p, q) = sMat->get(h, p, q);
                Sb(p + nc, q + nc) = sMat->get(h, p, q);

                Tb(p, q) = tMat->get(h, p, q);
                Tb(p + nc, q + nc) = tMat->get(h, p, q);

                Vb(p, q) = vMat->get(h, p, q);
                Vb(p + nc, q + nc) = vMat->get(h, p, q);

                Wb(p, q) = wMats[0]->get(h, p, q) + i * wMats[3]->get(h, p, q);
                Wb(p + nc, q + nc) = wMats[0]->get(h, p, q) - i * wMats[3]->get(h, p, q);
                Wb(p, q + nc) = wMats[2]->get(h, p, q) + i * wMats[1]->get(h, p, q);
                Wb(p + nc, q) = - wMats[2]->get(h, p, q) + i * wMats[1]->get(h, p, q);
            }
        }
    }
}

/*   Form the Dirac Hamiltonian
 *        | V       T       |
 *    d = |                 |
 *        | T  1/4c^2 W - T |
 */
void SOX2C1e::form_dirac_hamiltonian(SharedComplexMatrix Hdirac) {
    Hdirac->zero();

    auto tmp = std::make_shared<ComplexMatrix>("tmp", nsopi2c_.blocks());
    // recall that rel_pot_ has already been scaled by .25 / (pc_c_au * pc_c_au)
    (*tmp) = (*rel_pot_);
    (*tmp) -= (*kinetic_);

    for (int h = 0; h < nsopi2c_.n(); h++) {
        auto& Hblock = Hdirac->block(h);

        einsums::Range a{0, nsopi2c_[h]};
        einsums::Range b{nsopi2c_[h], 2*nsopi2c_[h]};

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
void SOX2C1e::form_orth(SharedComplexMatrix orth) {
    orth->zero();

    sOrth = sMat->clone();
    sOrth->power(-.5);

    tMat->scale(.5 / (pc_c_au * pc_c_au));
    tMat->power(-.5);

    for (int h = 0; h < nsopi_.n(); h++) {
        int dim = nsopi_[h];
        auto& orth_block = orth->block(h);
        for (int p = 0; p < dim; p++) {
            for (int q = 0; q < dim; q++) {
                orth_block(p, q) = sOrth->get(h, p, q);
                orth_block(p + dim, q + dim) = sOrth->get(h, p, q);
                orth_block(p + 2 * dim, q + 2 * dim) = tMat->get(h, p, q);
                orth_block(p + 3 * dim, q + 3 * dim) = tMat->get(h, p, q);
            }
        }
    }
}

void SOX2C1e::form_X(ComplexMatrix const& Hevec, SharedComplexMatrix X) {
    X->zero();
    for (int h = 0; h < nsopi_.n(); h++) {
        if (nsopi_[h] == 0) continue;
        einsums::Range small{2 * nsopi_[h], 4 * nsopi_[h]};
        auto pe = small;  // Positive energy solutions
        einsums::Range large{0, 2 * nsopi_[h]};

        einsums::TensorView C_large = Hevec[h](large, pe);
        einsums::TensorView C_small = Hevec[h](small, pe);

        // TODO: replace this with lapack call
        einsums::linear_algebra::invert(&C_large);
        gemm(C_small, C_large, &X->block(h));
    }
}

/* Stilde = S + \frac{1}{2c^2}X^\dagger T X */
void SOX2C1e::form_Stilde(ComplexMatrix const& X_HTX, SharedComplexMatrix Stilde) {
    const std::complex<double> C2{.5 / (pc_c_au * pc_c_au)};

    Stilde->zero();
    Stilde->operator+=(X_HTX);
    einsums::linear_algebra::scale(C2, Stilde.get());

    Stilde->operator+=(*overlap_);
}

/* Form R = S^{-1/2}[S^{-1/2} \tilde S S^{-1/2}]^{-1/2} S^{1/2} */
void SOX2C1e::form_R(ComplexMatrix const& Stilde, ComplexMatrix const& orth, SharedComplexMatrix R) {
    R->zero();
    SharedMatrix _shalf = sOrth->clone();
    _shalf->general_invert(); // Shalf == S^{1/2}

    auto orth2c = std::make_shared<ComplexMatrix>("orth top-left block", nsopi2c_.blocks());
    auto Shalf = std::make_shared<ComplexMatrix>("sMat^{1/2}", nsopi2c_.blocks());
    orth2c->zero();
    Shalf->zero();
    for (int h = 0; h < nsopi_.n(); h++) {
        auto& sblock = Shalf->block(h);
        for (int p = 0; p < nsopi_[h]; p++) {
            for (int q = 0; q < nsopi_[h]; q++) {
                Shalf->block(h)(p, q) = _shalf->get(h, p, q);
                Shalf->block(h)(p + nsopi_[h], q + nsopi_[h]) = _shalf->get(h, p, q);
            }
        }

        // Set 2c orth to top left of orth 4c
        einsums::Range a{0, nsopi2c_[h]};
        orth2c->block(h) = orth.block(h)(a, a);
    }

    auto tmp = std::make_shared<ComplexMatrix>("tmp", nsopi2c_.blocks());
    auto evec = std::make_shared<ComplexMatrix>("evec", nsopi2c_.blocks());
    auto eval = einsums::BlockTensor<double, 1>("eigenvalues", nsopi2c_.blocks());

    tmp->zero();
    evec->zero();

    // tmp <- [S^{-1/2} \tilde S S^{-1/2}]
    gemm(*orth2c, Stilde, evec);
    gemm(*evec, *orth2c, tmp);

    // Computing [S^{-1/2} \tilde S S^{-1/2}]^{-1/2}

    // TODO: switch to pow when einsums bug is fixed
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
        // Einsums needs block to be non-zero dim
        if (nsopi2c_[h] == 0) continue;
        auto& vecblock = evec->block(h);
        for (int i = 0; i < nsopi2c_[h]; i++) {
            einsums::linear_algebra::scale_column(i, std::pow(eval[h](i), -.5), &vecblock);
        }
    }
    gemm(*evec, *tmp, R);

    // S^{-1/2}...S^{1/2}
    gemm(*orth2c, *R, tmp);
    gemm(*tmp, *Shalf, R);
}

/* Takes Tx2c and Vx2c (two-component complex einsums BlockTensors) and converts
 * to two scalar matrices (real SharedMatrices T & V) and three Pauli components:
 *
 *     Tx2c + Vx2c ≈ T + V + iσ⋅(Hx, Hy, Hz)
 *
 * There may be a very small error in this conversion. Norm of this is printed.
 */
void SOX2C1e::form_pauli(ComplexMatrix& Tx2c, ComplexMatrix& Vx2c, SharedMatrix T, SharedMatrix V, SharedMatrix Hx,
                         SharedMatrix Hy, SharedMatrix Hz) {
    using namespace std::literals::complex_literals;
    double Verr = 0;
    double Terr = 0;
    double Xerr = 0;
    double Yerr = 0;
    double Zerr = 0;

    for (int h = 0; h < nsopi_.n(); h++) {
        int nc = nsopi_[h];
        auto& Tblock = Tx2c.block(h);
        auto& Vblock = Vx2c.block(h);

        einsums::Range a{0, nsopi_[h]};
        einsums::Range b{nsopi_[h], nsopi2c_[h]};

        for (int p = 0; p < nc; p++) {
            for (int q = 0; q < nc; q++) {
                std::complex<double>&& ts = .5 * (Tblock(a, a)(p, q) + Tblock(b, b)(p, q));
                Terr += ts.imag() * ts.imag();
                T->set(h, p, q, ts.real());

                std::complex<double>&& vs = .5 * (Vblock(a, a)(p, q) + Vblock(b, b)(p, q));
                Verr += vs.imag() * vs.imag();
                V->set(h, p, q, vs.real());

                std::complex<double>&& tz = -.5i * (Tblock(a, a)(p, q) - Tblock(b, b)(p, q));
                std::complex<double>&& vz = -.5i * (Vblock(a, a)(p, q) - Vblock(b, b)(p, q));
                Zerr += tz.imag() * tz.imag() + vz.imag() * vz.imag();
                Hz->set(h, p, q, tz.real() + vz.real());

                std::complex<double>&& ty = .5 * (Tblock(a, b)(p, q) - Tblock(a, b)(p, q));
                std::complex<double>&& vy = .5 * (Vblock(a, b)(p, q) - Vblock(a, b)(p, q));
                Yerr += ty.imag() * ty.imag() + vy.imag() * vy.imag();
                Hy->set(h, p, q, ty.real() + vy.real());

                std::complex<double>&& tx = -.5i * (Tblock(a, b)(p, q) + Tblock(b, a)(p, q));
                std::complex<double>&& vx = -.5i * (Vblock(a, b)(p, q) + Vblock(b, a)(p, q));
                Xerr += tx.imag() * tx.imag() + vx.imag() * vx.imag();
                Hx->set(h, p, q, tx.real() + vx.real());
            }
        }
    }

    outfile->Printf("  Frobenius norms of form_pauli\n");
    outfile->Printf("    T : %.13f\n", std::pow(Terr, .5));
    outfile->Printf("    V : %.13f\n", std::pow(Verr, .5));
    outfile->Printf("    Hx: %.13f\n", std::pow(Xerr, .5));
    outfile->Printf("    Hy: %.13f\n", std::pow(Yerr, .5));
    outfile->Printf("    Hz: %.13f\n", std::pow(Zerr, .5));
    if ((Terr > 1e-9) || (Verr > 1e-9) || (Xerr > 1e-9) || (Yerr > 1e-9) || (Zerr > 1e-9))
        outfile->Printf("  WARNING: Large error converting X2C ints to scalar and Pauli component form.");
}

/* Return the uncontracted to contracted basis transformation */
SharedMatrix SOX2C1e::get_projection() {
    // Integral factory for the BASIS/X2C_BASIS mixed basis
    auto integral_contracted = std::make_shared<IntegralFactory>(aoBasis_, deconBasis_, deconBasis_, deconBasis_);

    auto soBasis_contracted = std::make_shared<SOBasisSet>(aoBasis_, integral_contracted);

    nsopi_contracted_ = soBasis_contracted->dimension();

    auto soFactory_contracted = std::make_shared<MatrixFactory>();
    soFactory_contracted->init_with(nsopi_contracted_, nsopi_);

    // Form the overlap matrix in the BASIS/X2C_BASIS basis
    std::shared_ptr<OneBodySOInt> sOBI_cu(integral_contracted->so_overlap());
    SharedMatrix S_cu(soFactory_contracted->create_matrix("Overlap"));
    sOBI_cu->compute(S_cu);

    SharedMatrix S_inv = sMat->clone();
    S_inv->general_invert();

    auto D = std::make_shared<Matrix>("D", nsopi_, nsopi_contracted_);
    //  Form D = S_uu^{-1} S_uc.  Notice that we transpose S_cu
    D->gemm(false, true, 1.0, S_inv, S_cu, 0.0);
    return D;
}

/* Asserts that the X2C ints still reproduce the relevant 4c Dirac eigenvalues. */
bool SOX2C1e::test_hFW(einsums::BlockTensor<double, 1>& Heval, SharedMatrix S, SharedMatrix T, SharedMatrix V, SharedMatrix Hx,
              SharedMatrix Hy, SharedMatrix Hz) {
    using namespace std::literals::complex_literals;

    Dimension nsopi2c_contracted = nsopi_contracted_ + nsopi_contracted_;

    auto hFW = std::make_shared<ComplexMatrix>("Test hFW", nsopi2c_contracted.blocks());
    auto _orth = std::make_shared<Matrix>("orthogonalization", nsopi_contracted_, nsopi_contracted_);
    _orth = S->clone();
    _orth->power(-.5);
    hFW->zero();

    auto orth = std::make_shared<ComplexMatrix>("orthogonalization", nsopi2c_contracted.blocks());
    orth->zero();

    for (int h = 0; h < nsopi_contracted_.n(); h++) {
        int nc = nsopi_contracted_[h];

        auto& Ob = orth->block(h);
        auto& Hb = hFW->block(h);

        for (int p = 0; p < nc; p++) {
            for (int q = 0; q < nc; q++) {
                Hb(p, q)           += T->get(h, p, q) + V->get(h, p, q) + 1i * Hz->get(h, p, q);
                Hb(p, q + nc)      += Hx->get(h, p, q) + 1i * Hy->get(h, p, q);
                Hb(p + nc, q)      += Hx->get(h, p, q) - 1i * Hy->get(h, p, q);
                Hb(p + nc, q + nc) += T->get(h, p, q) + V->get(h, p, q) - 1i * Hz->get(h, p, q);

                Ob(p, q) += _orth->get(h, p, q);
                Ob(p + nc, q + nc) += _orth->get(h, p, q);
            }
        }
    }

    auto tmp = std::make_shared<ComplexMatrix>("tmp", nsopi2c_contracted.blocks());
    auto orth_H = std::make_shared<ComplexMatrix>("orth.conj().T", nsopi2c_contracted.blocks());
    conj_T(*orth, orth_H);

    gemm(*orth_H, *hFW, tmp);
    gemm(*tmp, *orth, hFW);

    auto test_eval = einsums::BlockTensor<double, 1>("hFW eigenvalues", nsopi2c_contracted.blocks());
    outfile->Printf("\n  Testing X2C Hamiltonian.\n");
    double sum = 0.0;
    bool print_eigenvalues = Process::environment.options.get_int("PRINT") > 1;
    for (int h = 0; h < nsopi_contracted_.n(); h++) {
        if (nsopi2c_contracted[h] == 0) continue;
        if (nsopi2c_contracted[h] > nsopi_[h])
            throw PSIEXCEPTION("SOX2C1e: RELATIVISTIC_BASIS is smaller than BASIS.");

        einsums::linear_algebra::heev<false>(&hFW->block(h), &test_eval[h]);

        outfile->Printf("    Irrep %d: Comparing (%d/%d) eigenvalues of Dirac Hamiltonian.\n", h+1, nsopi2c_contracted[h], nsopi_[h]);
        if (print_eigenvalues) outfile->Printf("      (eval)   |  H dirac  |   |   H X2C   |\n");
        for (int i = 0; i < nsopi2c_contracted[h]; i++) {
            sum += std::fabs(Heval[h](i + nsopi2c_[h]) - test_eval[h](i));

            if (print_eigenvalues) outfile->Printf("      (%4d)   %13.6f   %13.6f\n", i+1, Heval[h](i+nsopi2c_[h]), test_eval[h](i));
        }
    }

    outfile->Printf("\n  The 1-norm of |H_X2C - H_Dirac| is %.12f.\n", sum);
    return sum < 1.0e-6;
}

}  // namespace psi

#endif  // USING_Einsums
