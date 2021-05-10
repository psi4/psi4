/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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
FnoCC::FnoCC(std::string name, int naux, int nao, int nfrzc, int naocc, int nvir, 
             const SharedTensor2d& Corb, const SharedTensor2d& Fock, const SharedTensor2d& bQ, 
             double fno_cutoff, double fno_perct, int navir_fno,
             bool fno_by_perct, bool fno_by_user) {
    name_ = name;
    //ref_wfn_ = ref_wfn;
    //spin_ = spin;
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


FnoCC::FnoCC(std::string name,  
             int naux, int nao, int nfrzc, 
             int naoccA, int naoccB, int nvirA, int nvirB, 
             const SharedTensor2d& CorbA, const SharedTensor2d& CorbB,
             const SharedTensor2d& FockA, const SharedTensor2d& FockB,
             const SharedTensor2d& bQA, const SharedTensor2d& bQB, 
             double fno_cutoff, double fno_perct, int navir_fno,
             bool fno_by_perct, bool fno_by_user) 
{
    name_ = name;
    //ref_wfn_ = ref_wfn;
    naux_ = naux;
    naob_ = nao;
    naoccA_ = naoccA;
    naoccB_ = naoccB;
    nfrzc_ = nfrzc;
    noccA_ = nfrzc_ + naoccA_;
    noccB_ = nfrzc_ + naoccB_;
    nvirA_ = nvirA;
    nvirB_ = nvirB;
    norbA_ = noccA_ + nvirA_;
    norbB_ = noccB_ + nvirB_;
    navir_fnoA_ = navir_fno;
    navir_fnoA_ = navir_fno;
    fno_cutoff_ = fno_cutoff;
    fno_percentage_ = fno_perct;
    fno_by_perct_ = fno_by_perct;
    fno_by_user_ = fno_by_user;
    cutoff_ = 1.0E-10;

    // malloc
    CorbA_ = std::make_shared<Tensor2d>("FNO Alpha MO Coefficients", naob_, norbA_);
    CorbA_->copy(CorbA);
    CorbB_ = std::make_shared<Tensor2d>("FNO Beta MO Coefficients", naob_, norbB_);
    CorbB_->copy(CorbB);

    FockA_ = std::make_shared<Tensor2d>("FNO Alpha Fock matrix", norbA_, norbA_);
    FockA_->copy(FockA);
    FockB_ = std::make_shared<Tensor2d>("FNO Beta Fock matrix", norbB_, norbB_);
    FockB_->copy(FockB);

    bQA_ = std::make_shared<Tensor2d>("FNO DF_BASIS_CC B (Q|IA)", naux_, naoccA_, nvirA_);
    bQA_->copy(bQA);
    bQB_ = std::make_shared<Tensor2d>("FNO DF_BASIS_CC B (Q|ia)", naux_, naoccB_, nvirB_);
    bQB_->copy(bQB);

    GammaA_ = std::make_shared<Tensor2d>("OPDM VV Block", nvirA_, nvirA_);
    GammaB_ = std::make_shared<Tensor2d>("OPDM vv Block", nvirB_, nvirB_);

    TvvA_ = std::make_shared<Tensor2d>("FNO T matrix <V|V>", nvirA_, nvirA_);
    TvvB_ = std::make_shared<Tensor2d>("FNO T matrix <v|v>", nvirB_, nvirB_);
    TmatA_ = std::make_shared<Tensor2d>("FNO Alpha T matrix", norbA_, norbA_);
    TmatB_ = std::make_shared<Tensor2d>("FNO Beta T matrix", norbB_, norbB_);
    VmatA_ = std::make_shared<Tensor2d>("Alpha FNO Coefficients", naob_, norbA_);
    VmatB_ = std::make_shared<Tensor2d>("Beta FNO Coefficients", naob_, norbB_);

    // 1D
    diag_nA_ = std::make_shared<Tensor1d>("Alpha Diag G1", nvirB_);
    diag_nB_ = std::make_shared<Tensor1d>("Beta Diag G1", nvirB_);

    // fno
    //std::cout << "Before compute_fno\n" << std::endl;
    compute_fno_uhf();

}  //


FnoCC::~FnoCC() {
    // RHF
    if (Corb_) Corb_.reset();
    if (Fock_) Fock_.reset();
    if (Tvv_) Tvv_.reset();
    if (Tmat_) Tmat_.reset();
    if (Vmat_) Vmat_.reset();
    if(Gamma_) Gamma_.reset();
    if (bQ_) bQ_.reset();
    if(diag_n_) diag_n_.reset();

    // UHF
    if (CorbA_) CorbA_.reset();
    if (CorbB_) CorbB_.reset();
    if (FockA_) FockA_.reset();
    if (FockB_) FockB_.reset();
    if (TvvA_) TvvA_.reset();
    if (TvvB_) TvvB_.reset();
    if (TmatA_) TmatA_.reset();
    if (TmatB_) TmatB_.reset();
    if (VmatA_) VmatA_.reset();
    if (VmatB_) VmatB_.reset();
    if (GammaA_) GammaA_.reset();
    if (GammaB_) GammaB_.reset();
    if (bQA_) bQA_.reset();
    if (bQB_) bQB_.reset();
    if(diag_nA_) diag_nA_.reset();
    if(diag_nB_) diag_nB_.reset();


}  //

//===========================//
//=== RHF OPDM ==============//
//===========================//
void FnoCC::opdm_rhf() {
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

    //Gamma_->print();

}  //

//===========================//
//=== UHF OPDM ==============//
//===========================//
void FnoCC::opdm_uhf() {
    // tensors
    SharedTensor2d t2, L, K, M;

    // G_AB = 1/2\sum_{M,N,F} t_MN^FA t_MN^FB
    //t2AA_ump2_direct(t2);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA_, nvirA_, naoccA_, nvirA_);
    //tei_iajb_chem_directAA(L);
    L->gemm(true, false, bQA_, bQA_, 1.0, 0.0);
    M = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|AB>", naoccA_, naoccA_, nvirA_, nvirA_);
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ||AB>", naoccA_, naoccA_, nvirA_, nvirA_);
    //tei_pqrs_anti_symm_direct(K, M);
    K->sort(1243, M, 1.0, 0.0);
    K->scale(-1.0);
    K->add(M);
    M.reset();
    t2 = std::make_shared<Tensor2d>("T2_1 <IJ|AB>", naoccA_, naoccA_, nvirA_, nvirA_);
    t2->copy(K);
    t2->apply_denom(nfrzc_, noccA_, FockA_);
    Ecorr_ = 0.25 * t2->vector_dot(K);
    K.reset();

    GammaA_->contract442(4, 4, t2, t2, 0.5, 0.0);
    t2.reset();

    // G_AB += \sum_{M,n,f} t_Mn^Af t_Mn^Bf
    //t2AB_ump2_direct(t2);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|jb)", naoccA_, nvirA_, naoccB_, nvirB_);
    //tei_iajb_chem_directAB(L);
    L->gemm(true, false, bQA_, bQB_, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA_, naoccB_, nvirA_, nvirB_);
    K->sort(1324, L, 1.0, 0.0);
    L.reset();
    t2 = std::make_shared<Tensor2d>("T2_1 <Ij|Ab>", naoccA_, naoccB_, nvirA_, nvirB_);
    t2->copy(K);
    t2->apply_denom_os(nfrzc_, noccA_, noccB_, FockA_, FockB_);
    Ecorr_ += t2->vector_dot(K);
    K.reset();
    GammaA_->contract442(3, 3, t2, t2, 1.0, 1.0);

    // G_ab += \sum_{m,N,F} t_Mn^Fa t_Mn^Fb
    GammaB_->contract442(4, 4, t2, t2, 1.0, 0.0);
    t2.reset();

    // G_ab = -1/2\sum_{m,n,f} t_mn^fa t_mn^fb
    t2 = std::make_shared<Tensor2d>("T2_1 <ij|ab>", naoccB_, naoccB_, nvirB_, nvirB_);
    //t2BB_ump2_direct(t2);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ia|jb)", naoccB_, nvirB_, naoccB_, nvirB_);
    //tei_iajb_chem_directBB(L);
    L->gemm(true, false, bQB_, bQB_, 1.0, 0.0);
    M = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <ij|ab>", naoccB_, naoccB_, nvirB_, nvirB_);
    M->sort(1324, L, 1.0, 0.0);
    L.reset();
    K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <ij||ab>", naoccB_, naoccB_, nvirB_, nvirB_);
    //tei_pqrs_anti_symm_direct(K, M);
    K->sort(1243, M, 1.0, 0.0);
    K->scale(-1.0);
    K->add(M);
    M.reset();
    t2->copy(K);
    t2->apply_denom(nfrzc_, noccB_, FockB_);
    Ecorr_ += 0.25 * t2->vector_dot(K);
    K.reset();

    GammaB_->contract442(4, 4, t2, t2, 0.5, 1.0);
    t2.reset();

    //GammaA_->scale(2.0);
    //GammaB_->scale(2.0);
    //GammaA_->print();
    //GammaB_->print();

}//

//===========================//
//=== RHF FNO ===============//
//===========================//
void FnoCC::compute_fno() {
    SharedTensor2d K, T, U, X;
    SharedTensor1d diag;

    // OPDM
    opdm_rhf();

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

//===========================//
//=== UHF FNO ===============//
//===========================//
void FnoCC::compute_fno_uhf() {
    SharedTensor2d K, T, U, X;
    SharedTensor1d diag;

    // OPDM
    opdm_uhf();
    //std::cout << "UHF opdm is done\n" << std::endl;

    //=====================//
    // Alpha Part
    //=====================//
    // Diagonalize OPDM
    GammaA_->diagonalize(TvvA_, diag_nA_, cutoff_, false);
    GammaB_->diagonalize(TvvB_, diag_nB_, cutoff_, false);

    // determine frozen nos
    navirA_ = 0;
    navirB_ = 0;

    // if it is explicitly specified by user
    if (fno_by_user_) {
        navirA_ = navir_fno_;
        navirB_ = navir_fno_;
    }

    // if it is specified by percentage
    else if (fno_by_perct_) {
        double frac = fno_percentage_/100.0;
        navirA_ = (int)round(frac*nvirA_);
        navirB_ = (int)round(frac*nvirB_);
    }

    // if it is specified by OCC_TOLERANCE
    else {
        // Alpha 
        for (int a = 0; a < nvirA_; ++a) {
            double value = diag_nA_->get(a);
            if (std::fabs(diag_nA_->get(a)) >= fno_cutoff_) navirA_++;
        }

        // Beta
        for (int a = 0; a < nvirB_; ++a) {
            double value = diag_nB_->get(a);
            if (std::fabs(diag_nB_->get(a)) >= fno_cutoff_) navirB_++;
        }
    }
    nfrzvA_ = nvirA_ - navirA_;
    nfrzvB_ = nvirB_ - navirB_;

    // Equalize alpha-beta
    if (nfrzvB_ > nfrzvA_) {
        nfrzvB_ = nfrzvA_;
        navirB_ = nvirB_ - nfrzvB_;
    }
    else if (nfrzvB_ < nfrzvA_) {
        nfrzvA_ = nfrzvB_;
        navirA_ = nvirA_ - nfrzvA_;
    }
    nfrzv_ = nfrzvA_;


    // main if
    if (nfrzvA_ > 0) {
        // Semicanonic Fock
        SharedTensor2d UfvA = std::make_shared<Tensor2d>("UfvA", nfrzvA_, nfrzvA_);
        SharedTensor2d UvvA = std::make_shared<Tensor2d>("UvvA", navirA_, navirA_);
        SharedTensor2d FockfvA = std::make_shared<Tensor2d>("Fock <App|Bpp>", nfrzvA_, nfrzvA_);
        SharedTensor2d FockvvA = std::make_shared<Tensor2d>("Fock <Ap|Bp>", navirA_, navirA_);
        SharedTensor1d eigfvA = std::make_shared<Tensor1d>("Epsilon <App>", nfrzvA_);
        SharedTensor1d eigvvA = std::make_shared<Tensor1d>("Epsilon <Ap>", navirA_);

        // Form VV Block of old Fock matrix
        SharedTensor2d oldFvv_ = std::make_shared<Tensor2d>("Old Fock <V|V>", nvirA_, nvirA_);
        SharedTensor2d newFvv_ = std::make_shared<Tensor2d>("New Fock <V|V>", nvirA_, nvirA_);
#pragma omp parallel for
        for (int a = 0; a < nvirA_; ++a) {
            int aa = a + noccA_;
            for (int b = 0; b < nvirA_; ++b) {
                int bb = b + noccA_;
                oldFvv_->set(a, b, FockA_->get(aa, bb));
            }
        }
        newFvv_->transform(oldFvv_, TvvA_);
        oldFvv_.reset();

// set new VV block
#pragma omp parallel for
        for (int a = 0; a < nvirA_; ++a) {
            int aa = a + noccA_;
            for (int b = 0; b < nvirA_; ++b) {
                int bb = b + noccA_;
                FockA_->set(aa, bb, newFvv_->get(a, b));
            }
        }
        newFvv_.reset();
// Fock_->print();

// Fock_fv alpha spin case
#pragma omp parallel for
        for (int a = 0; a < nfrzvA_; ++a) {
            int aa = a + noccA_ + navirA_;
            for (int b = 0; b < nfrzvA_; ++b) {
                int bb = b + noccA_ + navirA_;
                FockfvA->set(a, b, FockA_->get(aa, bb));
            }
        }

// Fockvv alpha spin case
#pragma omp parallel for
        for (int a = 0; a < navirA_; ++a) {
            int aa = a + noccA_;
            for (int b = 0; b < navirA_; ++b) {
                int bb = b + noccA_;
                FockvvA->set(a, b, FockA_->get(aa, bb));
            }
        }

        // Diagonalize Fock
        FockfvA->diagonalize(UfvA, eigfvA, cutoff_);
        FockvvA->diagonalize(UvvA, eigvvA, cutoff_);

        // Build U
        SharedTensor2d Uorb_ = std::make_shared<Tensor2d>("FNO U orb", nvirA_, nvirA_);

// Ufv contribution alpha spin case
#pragma omp parallel for
        for (int a = 0; a < nfrzvA_; ++a) {
            int aa = a + navirA_;
            for (int b = 0; b < nfrzvA_; ++b) {
                int bb = b + navirA_;
                Uorb_->set(aa, bb, UfvA->get(a, b));
            }
        }

// Uvv contribution alpha spin case
#pragma omp parallel for
        for (int a = 0; a < navirA_; ++a) {
            for (int b = 0; b < navirA_; ++b) {
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
        SharedTensor2d Tvv_copy = std::make_shared<Tensor2d>("FNO T matrix copy", nvirA_, nvirA_);
        Tvv_copy->copy(TvvA_);
        TvvA_->gemm(false, false, Tvv_copy, Uorb_, 1.0, 0.0);
        Uorb_.reset();

    }  // end of main if

    // Form the overall T matrix
    for (int i = 0; i < noccA_; ++i) {
        TmatA_->set(i, i, 1.0);
    }

    for (int a = 0; a < nvirA_; ++a) {
        for (int b = 0; b < nvirA_; ++b) {
            TmatA_->set(a + noccA_, b + noccA_, TvvA_->get(a, b));
        }
    }

    // Form Vmatrix
    VmatA_->gemm(false, false, CorbA_, TmatA_, 1.0, 0.0);


    //=====================//
    // BETA Part
    //=====================//

    // main if
    if (nfrzvB_ > 0) {
        // Semicanonic Fock
        SharedTensor2d UfvB = std::make_shared<Tensor2d>("UfvB", nfrzvB_, nfrzvB_);
        SharedTensor2d UvvB = std::make_shared<Tensor2d>("UvvB", navirB_, navirB_);
        SharedTensor2d FockfvB = std::make_shared<Tensor2d>("Beta Fock <App|Bpp>", nfrzvB_, nfrzvB_);
        SharedTensor2d FockvvB = std::make_shared<Tensor2d>("Beta Fock <Ap|Bp>", navirB_, navirB_);
        SharedTensor1d eigfvB = std::make_shared<Tensor1d>("Beta Epsilon <App>", nfrzvB_);
        SharedTensor1d eigvvB = std::make_shared<Tensor1d>("Beta Epsilon <Ap>", navirB_);

        // Form VV Block of old Fock matrix
        SharedTensor2d oldFvv_ = std::make_shared<Tensor2d>("Old Fock <v|v>", nvirB_, nvirB_);
        SharedTensor2d newFvv_ = std::make_shared<Tensor2d>("New Fock <v|v>", nvirB_, nvirB_);
#pragma omp parallel for
        for (int a = 0; a < nvirB_; ++a) {
            int aa = a + noccB_;
            for (int b = 0; b < nvirB_; ++b) {
                int bb = b + noccB_;
                oldFvv_->set(a, b, FockB_->get(aa, bb));
            }
        }
        //std::cout << "fno ---> I am here\n" << std::endl;
        newFvv_->transform(oldFvv_, TvvB_);
        oldFvv_.reset();

// set new VV block
#pragma omp parallel for
        for (int a = 0; a < nvirB_; ++a) {
            int aa = a + noccB_;
            for (int b = 0; b < nvirB_; ++b) {
                int bb = b + noccB_;
                FockB_->set(aa, bb, newFvv_->get(a, b));
            }
        }
        newFvv_.reset();
// Fock_->print();

// Fock_fv alpha spin case
#pragma omp parallel for
        for (int a = 0; a < nfrzvB_; ++a) {
            int aa = a + noccB_ + navirB_;
            for (int b = 0; b < nfrzvB_; ++b) {
                int bb = b + noccB_ + navirB_;
                FockfvB->set(a, b, FockB_->get(aa, bb));
            }
        }

// Fockvv alpha spin case
#pragma omp parallel for
        for (int a = 0; a < navirB_; ++a) {
            int aa = a + noccB_;
            for (int b = 0; b < navirB_; ++b) {
                int bb = b + noccB_;
                FockvvB->set(a, b, FockB_->get(aa, bb));
            }
        }

        // Diagonalize Fock
        FockfvB->diagonalize(UfvB, eigfvB, cutoff_);
        FockvvB->diagonalize(UvvB, eigvvB, cutoff_);

        // Build U
        SharedTensor2d Uorb_ = std::make_shared<Tensor2d>("FNO U orb", nvirB_, nvirB_);

// Ufv contribution alpha spin case
#pragma omp parallel for
        for (int a = 0; a < nfrzvB_; ++a) {
            int aa = a + navirB_;
            for (int b = 0; b < nfrzvB_; ++b) {
                int bb = b + navirB_;
                Uorb_->set(aa, bb, UfvB->get(a, b));
            }
        }

// Uvv contribution alpha spin case
#pragma omp parallel for
        for (int a = 0; a < navirB_; ++a) {
            for (int b = 0; b < navirB_; ++b) {
                Uorb_->set(a, b, UvvB->get(a, b));
            }
        }

        UfvB.reset();
        UvvB.reset();
        FockfvB.reset();
        FockvvB.reset();
        eigfvB.reset();
        eigvvB.reset();

        // Get new Transformation matrix
        SharedTensor2d Tvv_copy = std::make_shared<Tensor2d>("Beta FNO T matrix copy", nvirB_, nvirB_);
        Tvv_copy->copy(TvvB_);
        TvvB_->gemm(false, false, Tvv_copy, Uorb_, 1.0, 0.0);
        Uorb_.reset();

    }  // end of main if

    //std::cout << "I am here\n" << std::endl;
    // Form the overall T matrix
    for (int i = 0; i < noccB_; ++i) {
        TmatB_->set(i, i, 1.0);
    }

    for (int a = 0; a < nvirB_; ++a) {
        for (int b = 0; b < nvirB_; ++b) {
            TmatB_->set(a + noccB_, b + noccB_, TvvB_->get(a, b));
        }
    }

    // Form Vmatrix
    VmatB_->gemm(false, false, CorbB_, TmatB_, 1.0, 0.0);

}  //


}  // namespace dfoccwave
}  // namespace psi
