/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libmints/matrix.h"

namespace psi {
namespace dfoccwave {

void DFOCC::ccd_t2_amps() {

    // RHF
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, I, T, Tnew, U, Tau, W, X, Y;
        SharedMatrix T2, RT2;

        // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
        // X(ia,jb) = \sum_{e} t_ij^ae F_be = \sum_{e} T(ia,je) F_be
        X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
        X->contract(false, true, naoccA * navirA * naoccA, navirA, navirA, t2, FabA, 1.0, 0.0);
        // X->cont424("IAJB", "IAJE", "BE", false, t2, FabA, 1.0, 0.0); // it works

        // t_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
        // X(ia,jb) = -\sum_{m} t_mj^ab F_mi = -\sum_{m} F(m,i) T(ma,jb)
        X->contract(true, false, naoccA, naoccA * navirA * navirA, naoccA, FijA, t2, -1.0, 1.0);
        // X->cont244("IAJB", "MI", "MAJB", false, FijA, t2, -1.0, 1.0); // it works
        X->symmetrize();

        // t_ij^ab <= <ij|ab>
        Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        Tnew->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

        // Contributions of X
        Tnew->axpy(X, 2.0);
        X.reset();

        // Write and close
        Tnew->write_symm(psio_, PSIF_DFOCC_AMPS);
        Tnew.reset();

        // WmnijT2
        ccd_WmnijT2();

        // WmbejT2
        ccd_WmbejT2();

        // WabefT2
        if (Wabef_type_ == "AUTO") {
            if (!do_ppl_hm)
                ccd_WabefT2();
            else
                ccd_WabefT2_high_mem();
        } else if (Wabef_type_ == "LOW_MEM")
            ccd_WabefT2();
        else if (Wabef_type_ == "HIGH_MEM")
            ccd_WabefT2_high_mem();

        // Denom
        Tnew = std::make_shared<Tensor2d>("New T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        Tnew->read_symm(psio_, PSIF_DFOCC_AMPS);
        Tnew->apply_denom_chem(nfrzc, noccA, FockA);

        // Reset T2
        rms_t2 = Tnew->rms(t2);

        // Error vector
        Tau = std::make_shared<Tensor2d>("RT2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        Tau->copy(Tnew);
        Tau->subtract(t2);
        t2->copy(Tnew);
        Tnew.reset();

        // DIIS
        RT2 = std::make_shared<Matrix>("RT2", naoccA * navirA, naoccA * navirA);
        Tau->to_matrix(RT2);
        Tau.reset();
        T2 = std::make_shared<Matrix>("T2", naoccA * navirA, naoccA * navirA);
        t2->to_matrix(T2);

        // add entry
        if (do_diis_ == 1 && orb_opt_ == "FALSE") ccsdDiisManager->add_entry(RT2.get(), T2.get());
        RT2.reset();

        // extrapolate
        if (do_diis_ == 1 && orb_opt_ == "FALSE") {
            if (ccsdDiisManager->subspace_size() >= cc_mindiis_) ccsdDiisManager->extrapolate(T2.get());
            t2->set2(T2);
        }
        T2.reset();

        // Form U(ia,jb) = 2*T(ia,jb) - T (ib,ja)
        U = std::make_shared<Tensor2d>("U2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        ccsd_u2_amps(U, t2);

        // Energy
        K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
        K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
        Ecorr = U->vector_dot(K);
        U.reset();
        K.reset();
        Eccd = Eref + Ecorr;
    }// if (reference_ == "RESTRICTED")

    // UHF
    else if (reference_ == "UNRESTRICTED") {
        ccd_W_MBEJAAAA();
        ccd_W_mbejBBBB();
        ccd_W_mBeJBABA();
        ccd_W_MbEjABAB();
        ccd_W_mBEjBAAB();
        ccd_W_MbeJABBA();
        //std::cout << "ccd_W_MBEJ is done \n";

        ccd_t2AA_amps();
        ccd_t2BB_amps();
        ccd_t2AB_amps();
        //std::cout << "ccd_t2 amps are done \n";
    }// else if (reference_ == "UNRESTRICTED")

}  // end ccd_t2_amps

void DFOCC::ccd_t2AA_amps()
{
    SharedTensor2d X, Y, K, L, T, T2, T2new, G, J, W, U, Tau, R;

    // t(IJ,AB) = <IJ||AB>
    J = std::make_shared<Tensor2d>("J (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <IJ||AB>", naoccA, naoccA, navirA, navirA);
    G->sort(1324, J, 1.0, 0.0);
    G->sort(1342, J, -1.0, 1.0);
    J.reset();

    // New T2AA
    T2new = std::make_shared<Tensor2d>("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2new->copy(G);
    G.reset();

    // t(IJ,AB) += P_(AB) \sum_(E) t(IJ,AE) * F_(B,E)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <IJ|AB>", naoccA, naoccA, navirA, navirA);
    X->contract(false, true, naoccA * naoccA * navirA, navirA, navirA, T2, FabA, 1.0, 0.0);
    T2.reset();
    T2new->add(X);
    T2new->sort(1243, X, -1.0, 1.0);
    X.reset();

    // t(IJ,AB) -= P_(IJ) \sum_(M) t(IM,AB) * F_(M,J)
    // t(IJ,AB) -= P_(IJ) \sum_(M) t(MJ,AB) * F_(M,I)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("T2 <MJ|AB>", naoccA, naoccA, navirA, navirA);
    Y->contract(true, false, naoccA, navirA * navirA * naoccA, naoccA, FijA, T2, -1.0, 0.0);
    T2.reset();
    T2new->sort(2134, Y, -1.0, 1.0);
    T2new->axpy(Y, 1.0);
    Y.reset();

    // t(IJ,AB) += P_(IJ)P_(AB) \sum_(M,E) T(IM,AE) * W_MBEJ
    // X(IA,JB) = \sum(ME) T(IA,ME) W(ME,JB)
    // T(IJ,AB) <= P_(IJ) * P_(AB) X(IA,JB)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    W = std::make_shared<Tensor2d>("W (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X1 (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->gemm(false, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();

    // t(IJ,AB) += P_(IJ)P_(AB) \sum_(m,e) T(Im,Ae) * W_mBeJ
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    U->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    X = std::make_shared<Tensor2d>("X2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    W = std::make_shared<Tensor2d>("W (me|JB)", naoccB, navirB, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(false, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();

    // Write and close
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new.reset();

    // WmnijT2
    ccd_W_MNIJT2AA();

    // WabefT2
    lccd_WabefT2AA();

    // Denom
    T2new = std::make_shared<Tensor2d>("New T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new->apply_denom(nfrzc, noccA, FockA);
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    //if (orb_opt_ == "FALSE") {
    if (do_diis_ == 1) {
        R = std::make_shared<Tensor2d>("RT2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        R->copy(T2new);
        R->subtract(T);
        R->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        R.reset();
    }
    T2new.reset();
    T.reset();
}// end ccd_t2AA_amps

void DFOCC::ccd_t2BB_amps()
{
    SharedTensor2d X, Y, K, L, T, T2, T2new, G, J, W, U, Tau, R;

    // t(ij,ab) = <ij||ab>
    J = std::make_shared<Tensor2d>("J (ia|jb)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <ij||ab>", naoccB, naoccB, navirB, navirB);
    G->sort(1324, J, 1.0, 0.0);
    G->sort(1342, J, -1.0, 1.0);
    J.reset();

    // New T2BB
    T2new = std::make_shared<Tensor2d>("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2new->copy(G);
    G.reset();

    // t(ij,ab) += P_(ab) \sum_(e) t(ij,ae) * F_(b,e)
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <ij|ab>", naoccB, naoccB, navirB, navirB);
    X->contract(false, true, naoccB * naoccB * navirB, navirB, navirB, T2, FabB, 1.0, 0.0);
    T2.reset();
    T2new->axpy(X, 1.0);
    T2new->sort(1243, X, -1.0, 1.0);
    X.reset();

    // t(ij,ab) -= P_(ij) \sum_(m) t(im,ab) * F_(m,j)
    // t(ij,ab) -= P_(ij) \sum_(m) t(mj,ab) * F_(m,i)
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("T2 <mj|ab>", naoccB, naoccB, navirB, navirB);
    Y->contract(true, false, naoccB, navirB * navirB * naoccB, naoccB, FijB, T2, -1.0, 0.0);
    T2.reset();
    T2new->sort(2134, Y, -1.0, 1.0);
    T2new->axpy(Y, 1.0);
    Y.reset();

    // t(ij,ab) += P_(ij) * P_(ab) \sum_(m,e) T(im,ae) * W_mbej
    // X(ia,jb) = \sum(m,e) T(ia,me) W(me,jb)
    // T(ij,ab) <= P_(ij) * P_(ab) X(ia,jb)
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 (ia|jb)", naoccB, navirB, naoccB, navirB);
    U->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    W = std::make_shared<Tensor2d>("W (me|jb)", naoccB, navirB, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X1 (ia|jb)", naoccB, navirB, naoccB, navirB);
    X->gemm(false, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();

    // t(ij,ab) += P_(ij) * P_(ab) \sum_(m,e) T(Mi,Ea) * W_MbEj
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    U->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    X = std::make_shared<Tensor2d>("X2 (ia|jb)", naoccB, navirB, naoccB, navirB);
    W = std::make_shared<Tensor2d>("W (ME|jb)", naoccA, navirA, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, U, W, 1.0, 0.0);
    U.reset();
    W.reset();
    //T2new->P_ijab(X);
    T2new->sort(1324, X, 1.0, 1.0);
    T2new->sort(3124, X, -1.0, 1.0);
    T2new->sort(1342, X, -1.0, 1.0);
    T2new->sort(3142, X, 1.0, 1.0);
    X.reset();

    // Write and close
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new.reset();

    // WmnijT2
    ccd_W_mnijT2BB();

    // WabefT2
    lccd_WabefT2BB();

    // Denom
    T2new = std::make_shared<Tensor2d>("New T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2new->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T2new->apply_denom(nfrzc, noccB, FockB);
    T2new->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    T = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    //if (orb_opt_ == "FALSE") {
    if (do_diis_ == 1) {
        R = std::make_shared<Tensor2d>("RT2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        R->copy(T2new);
        R->subtract(T);
        R->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        R.reset();
    }
    T2new.reset();
    T.reset();

}// End ccd_t2BB_amps

void DFOCC::ccd_t2AB_amps()
{
    SharedTensor2d X, Y, K, L, T, T2, T2new, G, J, W, U, Tau, R;

    // t(Ij,Ab) = <Ij|Ab>
    J = std::make_shared<Tensor2d>("J (IA|jb)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    T2new = std::make_shared<Tensor2d>("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2new->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Read T2
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);

    // T(Ij,Ab) += \sum_(e) T(Ij,Ae) F_(b,e)
    T2new->contract424(4, 2, T2, FabB, 1.0, 1.0);
    //T2new->contract(false, true, naoccA*naoccB*navirA, navirB, navirB, T2, FtabB, 1.0, 1.0);

    // T(Ij,Ab) += \sum_(E) T(Ij,Eb) Ft_(A,E)
    T2new->contract424(3, 2, T2, FabA, 1.0, 1.0);

    // T(Ij,Ab) -= \sum_(m) T(Im,Ab) Ft_(m,j)
    T2new->contract424(2, 1, T2, FijB, -1.0, 1.0);

    // T(Ij,Ab) -= \sum_(M) T(Mj,Ab) Ft_(M,I)
    T2new->contract424(1, 1, T2, FijA, -1.0, 1.0);
    T2.reset();

///WmbejT2AB Terms
    // T(ME,jb) <= T(Mj,Eb)
    // X(IA,jb) += \sum(ME) T(Mj,Eb) W(MA,EI) = \sum(ME) T(ME,jb) W(ME,IA)
    // T(Ij,Ab) <= X(IA,jb)
    W = std::make_shared<Tensor2d>("W (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    K->sort(1324, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X1 (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(true, false, W, K, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1324, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(IA,jb) = \sum(me) T(Im,Ae) W(mb,ej) = \sum(me) T(IA,me) W(me,jb)
    // T(Ij,Ab) <= X(IA,jb)
    W = std::make_shared<Tensor2d>("W (me|jb)", naoccB, navirB, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    K->sort(1324, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(false, false, K, W, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1324, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(IA,jb) = \sum(ME) T(IM,AE) W(Mb,Ej) = \sum(ME) T(IA,ME) W(ME,jb)
    // T(Ij,Ab) <= X(IA,jb)
    W = std::make_shared<Tensor2d>("W (ME|jb)", naoccA, navirA, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    K->sort(1324, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X3 (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(false, false, K, W, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1324, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(IA,jb) = \sum(me) T(jm,be) W(mA,eI) = \sum(me) T(me,jb) W(IA,me)
    // T(Ij,Ab) <= X(IA,jb)
    W = std::make_shared<Tensor2d>("W (me|JB)", naoccB, navirB, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (ia|jb)", naoccB, navirB, naoccB, navirB);
    K->sort(2413, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X4 (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(true, false, W, K, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1324, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(Ib,jA) = \sum(mE) T(Im,Eb) W(mA,Ej) = \sum(mE) T'(Ib,mE) W(mE,jA)
    // T(Ij,Ab) <= X(Ib,jA)
    X = std::make_shared<Tensor2d>("X5 (Ib|jA)", naoccA, navirB, naoccB, navirA);
    W = std::make_shared<Tensor2d>("W (mE|jB)", naoccB, navirA, naoccB, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (Ib|jA)", naoccA, navirB, naoccB, navirA);
    K->sort(1423, T2, 1.0, 0.0);
    T2.reset();

    X->gemm(false, false, K, W, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1342, X, 1.0, 1.0);
    X.reset();
    //T2new->write(psio_, PSIF_DFOCC_AMPS);
    //T2new.reset();

    // X(Ib,jA) = \sum(Me) T(Mj,Ae) W(Mb,eI) = \sum(me) T'(Me,jA) W(Me,Ib)
    // T(Ij,Ab) <= X(Ib,jA)
    W = std::make_shared<Tensor2d>("W (Me|Jb)", naoccA, navirB, naoccA, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("T2 (Ib|jA)", naoccA, navirB, naoccB, navirA);
    K->sort(1423, T2, 1.0, 0.0);
    T2.reset();

    X = std::make_shared<Tensor2d>("X6 (Ib|jA)", naoccA, navirB, naoccB, navirA);
    X->gemm(true, false, W, K, 1.0, 0.0);
    K.reset();
    W.reset();
    //T2new = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    //T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->sort(1342, X, 1.0, 1.0);
    X.reset();
    T2new->write(psio_, PSIF_DFOCC_AMPS);
    T2new.reset();

    // WmnijT2AB
    ccd_W_MnIjT2AB();
    //std::cout << "Wmnij \n";

    // WabefT2AB
    lccd_WabefT2AB();

    // Denom
    T2new = std::make_shared<Tensor2d>("New T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2new->read(psio_, PSIF_DFOCC_AMPS);
    T2new->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    T2new->write(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);

    //if (orb_opt_ == "FALSE" || mo_optimized == 1) {
    if (do_diis_ == 1) {
        R = std::make_shared<Tensor2d>("RT2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        R->copy(T2new);
        R->subtract(T);
        R->write(psio_, PSIF_DFOCC_AMPS);
        R.reset();
    }
    T2new.reset();
    T.reset();

    // DIIS
    //if (orb_opt_ == "FALSE" || mo_optimized == 1) {
    if (do_diis_ == 1) {
        SharedTensor2d RAA, TAA, RBB, TBB, RAB, TAB;
        // RAA
        RAA = std::make_shared<Tensor2d>("RT2 <IJ|AB>", ntri_anti_ijAA, ntri_anti_abAA);
        RAA->read(psio_, PSIF_DFOCC_AMPS);
        auto RT2AA = std::make_shared<Matrix>("RT2AA", ntri_anti_ijAA, ntri_anti_abAA);
        RAA->to_matrix(RT2AA);
        RAA.reset();
        // TAA
        TAA = std::make_shared<Tensor2d>("New T2 <IJ|AB>", ntri_anti_ijAA, ntri_anti_abAA);
        TAA->read(psio_, PSIF_DFOCC_AMPS);
        auto T2AA = std::make_shared<Matrix>("T2AA", ntri_anti_ijAA, ntri_anti_abAA);
        TAA->to_matrix(T2AA);

        // RBB
        RBB = std::make_shared<Tensor2d>("RT2 <ij|ab>", ntri_anti_ijBB, ntri_anti_abBB);
        RBB->read(psio_, PSIF_DFOCC_AMPS);
        auto RT2BB = std::make_shared<Matrix>("RT2BB", ntri_anti_ijBB, ntri_anti_abBB);
        RBB->to_matrix(RT2BB);
        RBB.reset();
        // TBB
        TBB = std::make_shared<Tensor2d>("New T2 <ij|ab>", ntri_anti_ijBB, ntri_anti_abBB);
        TBB->read(psio_, PSIF_DFOCC_AMPS);
        auto T2BB = std::make_shared<Matrix>("T2BB", ntri_anti_ijBB, ntri_anti_abBB);
        TBB->to_matrix(T2BB);

        // RAB
        RAB = std::make_shared<Tensor2d>("RT2 <Ij|Ab>", naoccA * naoccB, navirA * navirB);
        RAB->read(psio_, PSIF_DFOCC_AMPS);
        auto RT2AB = std::make_shared<Matrix>("RT2AB", naoccA * naoccB, navirA * navirB);
        RAB->to_matrix(RT2AB);
        RAB.reset();
        // TAB
        TAB = std::make_shared<Tensor2d>("New T2 <Ij|Ab>", naoccA * naoccB, navirA * navirB);
        TAB->read(psio_, PSIF_DFOCC_AMPS);
        auto T2AB = std::make_shared<Matrix>("T2AB", naoccA * naoccB, navirA * navirB);
        TAB->to_matrix(T2AB);

        // add entry
        //if (do_diis_ == 1)
            ccsdDiisManager->add_entry(RT2AA.get(), RT2BB.get(), RT2AB.get(), T2AA.get(), T2BB.get(), T2AB.get());
        RT2AA.reset();
        RT2BB.reset();
        RT2BB.reset();

        // extrapolate
        //if (do_diis_ == 1) {
            if (ccsdDiisManager->subspace_size() >= cc_mindiis_)
                ccsdDiisManager->extrapolate(T2AA.get(), T2BB.get(), T2AB.get());
            TAA->set2(T2AA);
            TBB->set2(T2BB);
            TAB->set2(T2AB);

        //}
        T2AA.reset();
        T2BB.reset();
        T2AB.reset();

        TAA->write(psio_, PSIF_DFOCC_AMPS);
        TBB->write(psio_, PSIF_DFOCC_AMPS);
        TAB->write(psio_, PSIF_DFOCC_AMPS);
        TAA.reset();
        TBB.reset();
        TAB.reset();
    }  // end if diis

    //=========================
    // Reset & Energy
    //=========================
    lccd_energy();
    Eccd = Elccd;
}// end ccsd_t2AB_amps


}  // namespace dfoccwave
}  // namespace psi
