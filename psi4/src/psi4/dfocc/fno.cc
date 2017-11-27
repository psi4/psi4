/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

// Latest revision on October 17, 2017.
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "fno.h"

namespace psi {
namespace dfoccwave {

/********************************************************************************************/
/************************** 1d array ********************************************************/
/********************************************************************************************/
FnoCC::FnoCC(std::string name, int naux, int nao, int nfrzc, int naocc, int nvir, const SharedTensor2d& Corb,
             const SharedTensor2d& Fock, const SharedTensor2d& bQ, double fno_cutoff, double fno_perct, int navir_fno,
             bool fno_by_perct, bool fno_by_user) {
    name_ = name;
    naux_ = naux;
    naob_ = nao;
    naocc_ = naocc;
    nfrzc_ = nfrzc;
    nocc_ = nfrzc_ + naocc_;
    nvir_ = nvir;
    norb_ = nocc_ + nvir_;
    navir_fno_ = navir_fno;
    fno_cutoff_ = fno_cutoff;
    fno_percentage_ = fno_perct;
    fno_by_perct_ = fno_by_perct;
    fno_by_user_ = fno_by_user;
    cutoff_ = 1.0E-10;

    // malloc
    Corb_ = std::make_shared<Tensor2d>("FNO Alpha MO Coefficients", naob_, norb_);
    Corb_->copy(Corb);
    Fock_ = std::make_shared<Tensor2d>("FNO Fock matrix", norb_, norb_);
    Fock_->copy(Fock);
    bQ_ = std::make_shared<Tensor2d>("FNO DF_BASIS_CC B (Q|IA)", naux_, naocc_, nvir_);
    bQ_->copy(bQ);
    Gamma_ = std::make_shared<Tensor2d>("OPDM VV Block", nvir_, nvir_);
    Tvv_ = std::make_shared<Tensor2d>("FNO T matrix <V|V>", nvir_, nvir_);
    Tmat_ = std::make_shared<Tensor2d>("FNO T matrix", norb_, norb_);
    Vmat_ = std::make_shared<Tensor2d>("Alpha FNO Coefficients", naob_, norb_);

    // 1D
    diag_n_ = std::make_shared<Tensor1d>("Diag G1", nvir_);

    // fno
    // std::cout << "Before compute_fno\n" << std::endl;
    compute_fno();

}  //

FnoCC::~FnoCC() {
    Corb_.reset();
    Fock_.reset();
    Tvv_.reset();
    Tmat_.reset();
    Vmat_.reset();
    Gamma_.reset();
    bQ_.reset();

    diag_n_.reset();
}  //

void FnoCC::compute_fno() {
    SharedTensor2d K, T, U, X;
    SharedTensor1d diag;

    // Compute T2
    // Build amplitudes in Mulliken order
    T = std::make_shared<Tensor2d>("FNO T2_1 (ia|jb)", naocc_, nvir_, naocc_, nvir_);
    K = std::make_shared<Tensor2d>("FNO (IA|JB)", naocc_, nvir_, naocc_, nvir_);
    K->gemm(true, false, bQ_, bQ_, 1.0, 0.0);
    T->copy(K);
    T->apply_denom_chem(nfrzc_, nocc_, Fock_);

    // form U(ia,jb)
    U = std::make_shared<Tensor2d>("FNO U2_1 (ia|jb)", naocc_, nvir_, naocc_, nvir_);
    U->sort(1432, T, 1.0, 0.0);
    U->scale(-1.0);
    U->axpy(T, 2.0);
    Ecorr_ = U->vector_dot(K);
    K.reset();

    // compute OPDM
    // G_ab = 2 \sum_{m,n,e} U'(ma,ne) T'(mb,ne)
    Gamma_->contract442(2, 2, U, T, 2.0, 0.0);
    T.reset();
    U.reset();

    // Diagonalize OPDM
    Gamma_->diagonalize(Tvv_, diag_n_, cutoff_, false);

    // determine frozen nos
    navir_ = 0;
    // if it is explicitly specified by user
    if (fno_by_user_) {
        navir_ = navir_fno_;
    }
    // if it is specified by percentage
    else if (fno_by_perct_) {
        double frac = fno_percentage_/100.0;
        navir_ = (int)round(frac*nvir_);
    }
    // if it is specified by OCC_TOLERANCE
    else {
        for (int a = 0; a < nvir_; ++a) {
            double value = diag_n_->get(a);
            if (std::fabs(diag_n_->get(a)) >= fno_cutoff_) navir_++;
        }
    }
    nfrzv_ = nvir_ - navir_;

    /*
    // check Unitary?
    X = std::make_shared<Tensor2d>("FNO Ttr*T matrix", norb_, norb_);
    X->gemm(true, false, Tmat_, Tmat_, 1.0, 0.0);
    X->print();
    X.reset();
    */

    // main if
    if (nfrzv_ > 0) {
        // Semicanonic Fock
        SharedTensor2d UfvA = std::make_shared<Tensor2d>("UfvA", nfrzv_, nfrzv_);
        SharedTensor2d UvvA = std::make_shared<Tensor2d>("UvvA", navir_, navir_);
        SharedTensor2d FockfvA = std::make_shared<Tensor2d>("Fock <App|Bpp>", nfrzv_, nfrzv_);
        SharedTensor2d FockvvA = std::make_shared<Tensor2d>("Fock <Ap|Bp>", navir_, navir_);
        SharedTensor1d eigfvA = std::make_shared<Tensor1d>("Epsilon <App>", nfrzv_);
        SharedTensor1d eigvvA = std::make_shared<Tensor1d>("Epsilon <Ap>", navir_);

        // Form VV Block of old Fock matrix
        SharedTensor2d oldFvv_ = std::make_shared<Tensor2d>("Old Fock <V|V>", nvir_, nvir_);
        SharedTensor2d newFvv_ = std::make_shared<Tensor2d>("New Fock <V|V>", nvir_, nvir_);
#pragma omp parallel for
        for (int a = 0; a < nvir_; ++a) {
            int aa = a + nocc_;
            for (int b = 0; b < nvir_; ++b) {
                int bb = b + nocc_;
                oldFvv_->set(a, b, Fock_->get(aa, bb));
            }
        }
        newFvv_->transform(oldFvv_, Tvv_);
        oldFvv_.reset();

// set new VV block
#pragma omp parallel for
        for (int a = 0; a < nvir_; ++a) {
            int aa = a + nocc_;
            for (int b = 0; b < nvir_; ++b) {
                int bb = b + nocc_;
                Fock_->set(aa, bb, newFvv_->get(a, b));
            }
        }
        newFvv_.reset();
// Fock_->print();

// Fock_fv alpha spin case
#pragma omp parallel for
        for (int a = 0; a < nfrzv_; ++a) {
            int aa = a + nocc_ + navir_;
            for (int b = 0; b < nfrzv_; ++b) {
                int bb = b + nocc_ + navir_;
                FockfvA->set(a, b, Fock_->get(aa, bb));
            }
        }

// Fockvv alpha spin case
#pragma omp parallel for
        for (int a = 0; a < navir_; ++a) {
            int aa = a + nocc_;
            for (int b = 0; b < navir_; ++b) {
                int bb = b + nocc_;
                FockvvA->set(a, b, Fock_->get(aa, bb));
            }
        }

        // Diagonalize Fock
        FockfvA->diagonalize(UfvA, eigfvA, cutoff_);
        FockvvA->diagonalize(UvvA, eigvvA, cutoff_);

        // Build U
        SharedTensor2d Uorb_ = std::make_shared<Tensor2d>("FNO U orb", nvir_, nvir_);

// Ufv contribution alpha spin case
#pragma omp parallel for
        for (int a = 0; a < nfrzv_; ++a) {
            int aa = a + navir_;
            for (int b = 0; b < nfrzv_; ++b) {
                int bb = b + navir_;
                Uorb_->set(aa, bb, UfvA->get(a, b));
            }
        }

// Uvv contribution alpha spin case
#pragma omp parallel for
        for (int a = 0; a < navir_; ++a) {
            for (int b = 0; b < navir_; ++b) {
                Uorb_->set(a, b, UvvA->get(a, b));
            }
        }

        UfvA.reset();
        UvvA.reset();
        FockfvA.reset();
        FockvvA.reset();
        eigfvA.reset();
        eigvvA.reset();

        // Get new Transformation matrix
        SharedTensor2d Tvv_copy = std::make_shared<Tensor2d>("FNO T matrix copy", nvir_, nvir_);
        Tvv_copy->copy(Tvv_);
        Tvv_->gemm(false, false, Tvv_copy, Uorb_, 1.0, 0.0);
        Uorb_.reset();

    }  // end of main if

    // Form the overall T matrix
    for (int i = 0; i < nocc_; ++i) {
        Tmat_->set(i, i, 1.0);
    }

    for (int a = 0; a < nvir_; ++a) {
        for (int b = 0; b < nvir_; ++b) {
            Tmat_->set(a + nocc_, b + nocc_, Tvv_->get(a, b));
        }
    }

    // Form Vmatrix
    Vmat_->gemm(false, false, Corb_, Tmat_, 1.0, 0.0);

}  //

}  // namespace dfoccwave
}  // namespace psi
