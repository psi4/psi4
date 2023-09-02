/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
#include "psi4/libmints/matrix.h"
#include "psi4/libdiis/diismanager.h"

namespace psi {
namespace dfoccwave {

void DFOCC::ccdl_l2_amps() {
    if (reference_ == "RESTRICTED") {
        // defs
        SharedTensor2d K, I, L, Lnew, T, U, Tau, W, X, Y, Z;

        // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
        // X(ia,jb) = \sum_{e} l_ij^ae F_eb = \sum_{e} L(ia,je) F_eb
        X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
        X->contract(false, false, naoccA * navirA * naoccA, navirA, navirA, l2, FabA, 1.0, 0.0);

        // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
        // X(ia,jb) = -\sum_{m} l_mj^ab F_im = -\sum_{m} F(i,m) L(ma,jb)
        X->contract(false, false, naoccA, naoccA * navirA * navirA, naoccA, FijA, l2, -1.0, 1.0);

        // l_ij^ab <= X(ia,jb) + X(jb,a) = 2Xt(ia,jb)
        // X(ia,jb) = \sum_{Q} (G_ai^Q - G_ia^Q) b_jb^Q
        U = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
        T->swap_3index_col(U);
        U.reset();
        U = std::make_shared<Tensor2d>("G (Q|IA)", nQ, naoccA, navirA);
        U->read(psio_, PSIF_DFOCC_AMPS);
        T->axpy(U, -1.0);
        U.reset();
        X->gemm(true, false, T, bQiaA, 1.0, 1.0);
        T.reset();
        X->symmetrize();

        // l_ij^ab <= <ij|ab>
        Lnew = std::make_shared<Tensor2d>("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        Lnew->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);

        // Contributions of X
        Lnew->axpy(X, 2.0);
        X.reset();

        // Write and close
        Lnew->write_symm(psio_, PSIF_DFOCC_AMPS);
        Lnew.reset();

        // VmnijL2
        ccdl_VmnijL2();

        // WijmnL2
        ccdl_WijmnL2();

        // WmbejTL2
        ccdl_WmbejL2();

        // WabefL2
        if (Wabef_type_ == "AUTO") {
            if (!do_ppl_hm)
                ccdl_WabefL2();
            else
                ccsdl_WabefL2_high_mem();
        } else if (Wabef_type_ == "LOW_MEM")
            ccdl_WabefL2();
        else if (Wabef_type_ == "HIGH_MEM")
            ccsdl_WabefL2_high_mem();

        // Denom
        Lnew = std::make_shared<Tensor2d>("New L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        Lnew->read_symm(psio_, PSIF_DFOCC_AMPS);
        Lnew->apply_denom_chem(nfrzc, noccA, FockA);

        // Reset L2
        rms_l2 = Lnew->rms(l2);

        // Error vector
        Tau = std::make_shared<Tensor2d>("RL2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        Tau->copy(Lnew);
        Tau->subtract(l2);
        l2->copy(Lnew);
        Lnew.reset();

        // DIIS
        auto RL2 = std::make_shared<Matrix>("RL2", naoccA * navirA, naoccA * navirA);
        Tau->to_matrix(RL2);
        Tau.reset();
        auto L2 = std::make_shared<Matrix>("L2", naoccA * navirA, naoccA * navirA);
        l2->to_matrix(L2);

        // add entry
        if (do_diis_ == 1 && orb_opt_ == "FALSE") ccsdlDiisManager->add_entry(RL2.get(), L2.get());
        RL2.reset();

        // extrapolate
        if (do_diis_ == 1 && orb_opt_ == "FALSE") {
            if (ccsdlDiisManager->subspace_size() >= cc_mindiis_) ccsdlDiisManager->extrapolate(L2.get());
            l2->set2(L2);
        }
        L2.reset();

        // Energy
        if (orb_opt_ == "FALSE") {
            U = std::make_shared<Tensor2d>("2*L(ia,jb) - L(ib,ja)", naoccA, navirA, naoccA, navirA);
            U->sort(1432, l2, 1.0, 0.0);
            U->scale(-1.0);
            U->axpy(l2, 2.0);
            K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
            K->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
            EcorrL = U->vector_dot(K);
            U.reset();
            K.reset();
            EccdL = Escf + EcorrL;
        }
    }// end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        ccdl_l2AA_amps();
        ccdl_l2BB_amps();
        ccdl_l2AB_amps();
    }// end else if (reference_ == "UNRESTRICTED")

}  // end ccdl_l2_amps

//======================================================================
//    UHF L2AA
//======================================================================
void DFOCC::ccdl_l2AA_amps()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L, U, V, Tau, Lnew, G, R, T2;
    // l_IJ^AB = <IJ||AB> + P_(AB) \sum_{E} l_IJ^AE F_EB - P_(IJ) \sum_{M} l_IM^AB F_JM + 1/2 \sum_{MN} l_MN^AB W_IJMN + 1/2 \sum_{EF} l_IJ^EF W_EFAB
    //         + \sum_{MN} V_MNIJ <MN|AB> +   P_(IJ) P_(AB) \sum_{ME} l_IM^AE W_JEBM
    //         +  P_(IJ) P_(AB) \sum_{me} l_Im^Ae W_JeBm + P_(IJ) P_(AB) \sum_{Q} (G_AI^Q - G_IA^Q) b_JB^Q                           (117)

    // l_IJ^AB += <IJ||AB>    (1)
    J = std::make_shared<Tensor2d>("J (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    Lnew = std::make_shared<Tensor2d>("New L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Lnew->sort(1324, J, 1.0, 0.0);
    Lnew->sort(1342, J, -1.0, 1.0);
    J.reset();
    // l_IJ^AB += P_(AB) \sum_{E} l_IJ^AE F_EB   (2)
    L = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <IJ|AB>", naoccA, naoccA, navirA, navirA);
    X->contract(false, false, naoccA * naoccA * navirA, navirA, navirA, L, FabA, 1.0, 0.0);
    Lnew->axpy(X, 1.0);
    Lnew->sort(1243, X, -1.0, 1.0);
    X.reset();
    // l_IJ^AB -= P_(IJ) \sum_{M} l_IM^AB F_JM   (3)
    Y = std::make_shared<Tensor2d>("Y <AB|IJ>", navirA, navirA, naoccA, naoccA);
    Y->trans(L);
    X = std::make_shared<Tensor2d>("X <AB|IJ>", navirA, navirA, naoccA, naoccA);
    X->contract(false, true, navirA * navirA * naoccA, naoccA, naoccA, Y, FijA, 1.0, 0.0);
    Y.reset();
    Lnew->sort(3412, X, -1.0, 1.0);
    Lnew->sort(4312, X, 1.0, 1.0);
    X.reset();
    // l_IJ^AB += 1/2 \sum_{MN} l_MN^AB W_IJMN    (4)
    W = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    Lnew->gemm(false, false, W, L, 0.5, 1.0);
    W.reset();
    L.reset();
    Lnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    // l_IJ^AB += \sum_{EF} l_IJ^EF W_EFAB    (5)
    cc_WabefT2AA("L");

    // l_IJ^AB += \sum_{MN} V_MNIJ <MN|AB>        (6)
    J = std::make_shared<Tensor2d>("J (MA|NB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <MN|AB>", naoccA, naoccA, navirA, navirA);
    G->sort(1324, J, 1.0, 0.0);
    J.reset();
    V = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Lnew = std::make_shared<Tensor2d>("New L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Lnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->gemm(true, false, V, G, 1.0, 1.0);
    V.reset();
    G.reset();

    // l_IJ^AB += P_(IJ) P_(AB) \sum_{ME} l_IM^AE W_JEBM    (9)
    Y = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Y->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 (ME|IA)", naoccA, navirA, naoccA, navirA);
    L->sort(2413, Y, 1.0, 0.0);
    Y.reset();
    W = std::make_shared<Tensor2d>("WL (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (JB|IA)", naoccA, navirA, naoccA, navirA);
    X->gemm(false, false, W, L, 1.0, 0.0);
    L.reset();
    W.reset();
    //Lnew->P_ijab(X);
    Lnew->sort(3142, X, 1.0, 1.0);
    Lnew->sort(1342, X, -1.0, 1.0);
    Lnew->sort(3124, X, -1.0, 1.0);
    Lnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // l_IJ^AB += P_(IJ) P_(AB) \sum_{me} l_Im^Ae W_JeBm    (10)
    Y = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 (me|IA)", naoccB, navirB, naoccA, navirA);
    L->sort(2413, Y, 1.0, 0.0);
    Y.reset();
    W = std::make_shared<Tensor2d>("WL (ME|jb)", naoccA, navirA, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (JB|IA)", naoccA, navirA, naoccA, navirA);
    X->gemm(false, false, W, L, 1.0, 0.0);
    L.reset();
    W.reset();
    //Lnew->P_ijab(X);
    Lnew->sort(3142, X, 1.0, 1.0);
    Lnew->sort(1342, X, -1.0, 1.0);
    Lnew->sort(3124, X, -1.0, 1.0);
    Lnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // l_IJ^AB += P_(IJ) P_(AB) \sum_{Q} (G_AI^Q - G_IA^Q) b_JB^Q      (11)
    U = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("G (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    X = std::make_shared<Tensor2d>("X (IA|JB)", naoccA, navirA, naoccA, navirA);
    X->gemm(true, false, T, bQiaA, 1.0, 0.0);
    T.reset();
    //Lnew->P_ijab(X);
    Lnew->sort(1324, X, 1.0, 1.0);
    Lnew->sort(3124, X, -1.0, 1.0);
    Lnew->sort(1342, X, -1.0, 1.0);
    Lnew->sort(3142, X, 1.0, 1.0);
    X.reset();

    // Denom
    Lnew->apply_denom(nfrzc, noccA, FockA);
    // Write and close
    Lnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    //Lnew->print();

    // Error vector
    L = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    R = std::make_shared<Tensor2d>("RL2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    R->copy(Lnew);
    Lnew.reset();
    R->subtract(L);
    L.reset();
    R->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
} // end ccdl_l2AA_amps

//======================================================================
//    UHF L2BB
//======================================================================
void DFOCC::ccdl_l2BB_amps()
{
    SharedTensor2d G, J, W, I, K, X, Y, T, Z, L, U, V, Tau, Lnew, R, T2;
    // l_ij^ab = <ij||ab> + P_(ab) \sum_{e} l_ij^ae F_eb - P_(ij) \sum_{m} l_im^ab F_jm + 1/2 \sum_{mn} l_mn^ab W_ijmn + 1/2 \sum_{ef} l_ij^ef W_efab
    //         + \sum_{mn} V_mnij <mn|ab> + P_(ij) P_(ab) \sum_{me} l_im^ae W_jebm
    //         + P_(ij) P_(ab) \sum_{ME} l_Mi^Ea W_jEbM + P_(ij) P_(ab) \sum_{Q} (G_ai^Q - G_ia^Q) b_jb^Q                             (118)

    // l_ij^ab +=  <ij||ab>    (1)
    J = std::make_shared<Tensor2d>("J (ia|jb)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    Lnew = std::make_shared<Tensor2d>("New L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    Lnew->sort(1324, J, 1.0, 0.0);
    Lnew->sort(1342, J, -1.0, 1.0);
    J.reset();
    // l_ij^ab +=  P_(ab) \sumi_{e} l_ij^ae F_eb   (2)
    L = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <ij|ab>", naoccB, naoccB, navirB, navirB);
    X->contract(false, false, naoccB * naoccB * navirB, navirB, navirB, L, FabB, 1.0, 0.0);
    Lnew->axpy(X, 1.0);
    Lnew->sort(1243, X, -1.0, 1.0);
    X.reset();
    // l_ij^ab -=  P_(ij) \sum_{m} l_im^ab F_jm   (3)
    Y = std::make_shared<Tensor2d>("Y <ab|ij>", navirB, navirB, naoccB, naoccB);
    Y->trans(L);
    X = std::make_shared<Tensor2d>("X <ab|ij>", navirB, navirB, naoccB, naoccB);
    X->contract(false, true, navirB * navirB * naoccB, naoccB, naoccB, Y, FijB, 1.0, 0.0);
    Y.reset();
    Lnew->sort(3412, X, -1.0, 1.0);
    Lnew->sort(4312, X, 1.0, 1.0);
    X.reset();
    // l_ij^ab +=  1/2 \sum_{mn} l_mn^ab W_ijmn    (4)
    W = std::make_shared<Tensor2d>("W <mn|ij>", naoccB, naoccB, naoccB, naoccB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    Lnew->gemm(false, false, W, L, 0.5, 1.0);
    W.reset();
    L.reset();
    Lnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    // l_ij^ab +=  \sum_{ef} l_ij^ef W_efab    (5)
    cc_WabefT2BB("L");

    // l_ij^ab +=  \sum_{mn} V_mnij <mn|ab>        (6)
    J = std::make_shared<Tensor2d>("J (ma|nb)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <mn|ab>", naoccB, naoccB, navirB, navirB);
    G->sort(1324, J, 1.0, 0.0);
    J.reset();
    V = std::make_shared<Tensor2d>("V <ij|kl>", naoccB, naoccB, naoccB, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Lnew = std::make_shared<Tensor2d>("New L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    Lnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Lnew->gemm(true, false, V, G, 1.0, 1.0);
    V.reset();
    G.reset();

    // l_ij^ab +=  P_(ij) P_(ab) \sum_{me} l_im^ae W_jebm   (9)
    Y = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    Y->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 (me|ia)", naoccB, navirB, naoccB, navirB);
    L->sort(2413, Y, 1.0, 0.0);
    Y.reset();
    W = std::make_shared<Tensor2d>("WL (me|jb)", naoccB, navirB, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (jb|ia)", naoccB, navirB, naoccB, navirB);
    X->gemm(false, false, W, L, 1.0, 0.0);
    L.reset();
    W.reset();
    //Lnew->P_ijab(X);
    Lnew->sort(3142, X, 1.0, 1.0);
    Lnew->sort(1342, X, -1.0, 1.0);
    Lnew->sort(3124, X, -1.0, 1.0);
    Lnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // l_ij^ab +=  P_(ij) P_(ab) \sum_{ME} l_Mi^Ea W_jEbM   (10)
    Y = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 (ME|ia)", naoccA, navirA, naoccB, navirB);
    L->sort(1324, Y, 1.0, 0.0);
    Y.reset();
    W = std::make_shared<Tensor2d>("WL (me|JB)", naoccB, navirB, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (jb|ia)", naoccB, navirB, naoccB, navirB);
    X->gemm(false, false, W, L, 1.0, 0.0);
    L.reset();
    W.reset();
    //Lnew->P_ijab(X);
    Lnew->sort(3142, X, 1.0, 1.0);
    Lnew->sort(1342, X, -1.0, 1.0);
    Lnew->sort(3124, X, -1.0, 1.0);
    Lnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // l_ij^ab +=  P_(ij) P_(ab) \sum_{Q} (G_ai^Q - G_ia^Q) b_jb^Q   (11)
    U = std::make_shared<Tensor2d>("G (Q|ai)", nQ, navirB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|ia)", nQ, naoccB, navirB);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("G (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    X = std::make_shared<Tensor2d>("X (ia|jb)", naoccB, navirB, naoccB, navirB);
    X->gemm(true, false, T, bQiaB, 1.0, 0.0);
    T.reset();
    //Lnew->P_ijab(X);
    Lnew->sort(1324, X, 1.0, 1.0);
    Lnew->sort(3124, X, -1.0, 1.0);
    Lnew->sort(1342, X, -1.0, 1.0);
    Lnew->sort(3142, X, 1.0, 1.0);
    X.reset();

    // Denom
    Lnew->apply_denom(nfrzc, noccB, FockB);
    // Write and close
    Lnew->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    //Lnew->print();

    // Error vector
    L = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    R = std::make_shared<Tensor2d>("RL2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    R->copy(Lnew);
    Lnew.reset();
    R->subtract(L);
    L.reset();
    R->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
} // end ccsdl_l2BB_amps

//======================================================================
//    UHF L2AB
//======================================================================
void DFOCC::ccdl_l2AB_amps()
{
    SharedTensor2d G, J, W, I, K, X, Y, T, Z, L, U, V, Tau, Lnew, R, T2;
    // l_Ij^Ab = <Ij|Ab> + \sum_{e} l_Ij^Ae Ft_eb + \sum_{E} l_Ij^Eb F_EA - \sum_{m} l_Im^Ab F_jm - \sum_{M} l_Mj^Ab Ft_IM + \sum_{Mn} l_Mn^Ab W_IjMn
    //         + \sum_{Ef} l_Ij^Ef W_EfAb + \sum_{Mn} V_MnIj <Mn|Ab>
    //         + \sum_{ME} l_IM^AE W_jEbM + \sum_{me} l_Im^Ae W_jebm + \sum_{Me} l_Mj^Ae W_IebM + \sum_{mE} l_Im^Eb W_jEAm + \sum_{me} l_jm^be W_IeAm
    //         + \sum_{ME} l_Mj^Eb W_IEAM + \sum_{Q} (G_AI^Q - G_IA^Q) b_jb^Q + \sum_{Q} (G_bj^Q - G_jb^Q) b_IA^Q   (119)

    // l_Ij^Ab += <Ij|Ab>    (1)
    J = std::make_shared<Tensor2d>("J (IA|jb)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    Lnew = std::make_shared<Tensor2d>("New L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Lnew->sort(1324, J, 1.0, 0.0);
    J.reset();

    // l_Ij^Ab += \sum_{e} l_Ij^Ae F_eb   (2)
    L = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    Lnew->contract(false, false, naoccA * naoccB * navirA, navirB, navirB, L, FabB, 1.0, 1.0);

    // l_Ij^Ab += \sum_{E} l_Ij^Eb F_EA   (3)
    Y = std::make_shared<Tensor2d>("L <Eb|Ij>", navirA, navirB, naoccA, naoccB);
    Y->trans(L);
    L.reset();
    X = std::make_shared<Tensor2d>("X <Ab|Ij>", navirA, navirB, naoccA, naoccB);
    X->contract(true, false, navirA, navirB * naoccA * naoccB, navirA, FabA, Y, 1.0, 0.0);

    // l_Ij^Ab -= \sum_{m} l_Im^Ab F_jm   (4)
    X->contract(false, true, navirA * navirB * naoccA, naoccB, naoccB, Y, FijB, -1.0, 1.0);
    Y.reset();
    Lnew->sort(3412, X, 1.0, 1.0);
    X.reset();

    // l_Ij^Ab -= \sum_{M} l_Mj^Ab F_IM   (5)
    L = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    Lnew->contract(false, false, naoccA, naoccB * navirA * navirB, naoccA, FijA, L, -1.0, 1.0);
    L.reset();

    // l_Ij^Ab += \sum_{Mn} l_Mn^Ab W_IjMn     (6)
    X = std::make_shared<Tensor2d>("W <Mn|Ij>", naoccA, naoccB, naoccA, naoccB);
    X->read(psio_, PSIF_DFOCC_AMPS);
    //W = std::make_shared<Tensor2d>("W <Ij|Mn>", naoccA, naoccB, naoccA, naoccB);
    //W->trans(X);
    //X.reset();
    L = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    Lnew->gemm(false, false, X, L, 1.0, 1.0);
    L.reset();
    X.reset();
    Lnew->write(psio_, PSIF_DFOCC_AMPS);
    Lnew.reset();

    // l_Ij^Ab += \sum_{Ef} l_Ij^Ef W_EfAb     (7)
    cc_WabefT2AB("L");

    // l_Ij^Ab += \sum_{Mn} V_MnIj <Mn|Ab>     (8)
    J = std::make_shared<Tensor2d>("J (MA|nb)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G <Mn|Ab>", naoccA, naoccB, navirA, navirB);
    G->sort(1324, J, 1.0, 0.0);
    J.reset();
    V = std::make_shared<Tensor2d>("V <Ij|Kl>", naoccA, naoccB, naoccA, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Lnew = std::make_shared<Tensor2d>("New L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Lnew->read(psio_, PSIF_DFOCC_AMPS);
    Lnew->gemm(true, false, V, G, 1.0, 1.0);
    V.reset();
    G.reset();

    // l_Ij^Ab += \sum_{ME} l_IM^AE W_jEbM    (13)
    U = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    L->sort(1324, U, 1.0, 0.0);
    U.reset();
    W = std::make_shared<Tensor2d>("WL (me|JB)", naoccB, navirB, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(false, true, L, W, 1.0, 0.0); // X acik kalcak asagida kullaniyorum
    L.reset();
    W.reset();

    // l_Ij^Ab += \sum_{me} l_Im^Ae W_jebm    (14)
    U = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 (IA|jb)", naoccA, navirA, naoccB, navirB);
    L->sort(1324, U, 1.0, 0.0);
    // U acik kalcak asagida kullaniyorum
    W = std::make_shared<Tensor2d>("WL (me|jb)", naoccB, navirB, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(false, true, L, W, 1.0, 1.0);
    L.reset();
    W.reset();
    Lnew->sort(1324, X, 1.0, 1.0);
    //UB start
    /*
    Y = std::make_shared<Tensor2d>("Y <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Y->sort(1324, X, 1.0, 0.0);
    Lnew->axpy(Y,1.0);
    Y.reset();
    */
    //UB end
    X.reset();

    // l_Ij^Ab += \sum_{Me} l_Mj^Ae W_IebM    (15)
    L = std::make_shared<Tensor2d>("L2 <Me|jA>", naoccA, navirB, naoccB, navirA);
    L->sort(1423, U, 1.0, 0.0);
    W = std::make_shared<Tensor2d>("WL (Me|Jb)", naoccA, navirB, naoccA, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (Ib|jA)", naoccA, navirB, naoccB, navirA);
    X->gemm(false, false, W, L, 1.0, 0.0);
    L.reset();
    W.reset();

    // l_Ij^Ab += \sum_{mE} l_Im^Eb W_jEAm    (16)
    L = std::make_shared<Tensor2d>("L2p <Ib|mE>", naoccA, navirB, naoccB, navirA);     ///////////////////// hata vardı yeni sort ile düzelttim
    L->sort(1423, U, 1.0, 0.0);
    W = std::make_shared<Tensor2d>("WL (mE|jB)", naoccB, navirA, naoccB, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(false, true, L, W, 1.0, 1.0);
    L.reset();
    W.reset();
    Lnew->sort(1342, X, 1.0, 1.0);
    X.reset();

    // l_Ij^Ab += \sum_{ME} l_Mj^Eb W_IEAM    (18)
    W = std::make_shared<Tensor2d>("WL (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 (ME|jb)", naoccA, navirA, naoccB, navirB);
    L->sort(1324, U, 1.0, 0.0);
    U.reset();
    X = std::make_shared<Tensor2d>("X (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(false, false, W, L, 1.0, 0.0);
    L.reset();
    W.reset();

    // l_Ij^Ab += \sum_{me} l_jm^be W_IeAm    (17)
    W = std::make_shared<Tensor2d>("WL (ME|jb)", naoccA, navirA, naoccB, navirB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    U->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 (me|jb)", naoccB, navirB, naoccB, navirB);
    L->sort(1324, U, 1.0, 0.0);
    U.reset();
    X->gemm(false, true, W, L, 1.0, 1.0);
    L.reset();
    W.reset();
    Lnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // l_Ij^Ab += \sum_{Q} (G_AI^Q - G_IA^Q) b_jb^Q    (19)
    U = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("G (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    X = std::make_shared<Tensor2d>("X (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(true, false, T, bQiaB, 1.0, 0.0);
    T.reset();
    Lnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // l_Ij^Ab += \sum_{Q} (G_bj^Q - G_jb^Q) b_IA^Q    (20)
    U = std::make_shared<Tensor2d>("G (Q|ai)", nQ, navirB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|ia)", nQ, naoccB, navirB);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("G (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    X = std::make_shared<Tensor2d>("X (IA|jb)", naoccA, navirA, naoccB, navirB);
    X->gemm(true, false, bQiaA, T, 1.0, 0.0);
    T.reset();
    Lnew->sort(1324, X, 1.0, 1.0);
    X.reset();

    // Denom
    Lnew->apply_denom_os(nfrzc, noccA, noccB, FockA, FockB);
    // Write and close
    Lnew->write(psio_, PSIF_DFOCC_AMPS);

    // Error vector
    L = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    R = std::make_shared<Tensor2d>("RL2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    R->copy(Lnew);
    Lnew.reset();
    R->subtract(L);
    L.reset();
    R->write(psio_, PSIF_DFOCC_AMPS);
    R.reset();


    // DIIS
    //if (orb_opt_ == "FALSE" && do_diis_ == 1) {
    if (do_diis_ == 1) {
        SharedTensor2d RAA, LAA, RBB, LBB, RAB, LAB;

            // RAA
            RAA = std::make_shared<Tensor2d>("RL2 <IJ|AB>", ntri_anti_ijAA, ntri_anti_abAA);
            RAA->read(psio_, PSIF_DFOCC_AMPS);
            auto RL2AA = std::make_shared<Matrix>("RL2AA", ntri_anti_ijAA, ntri_anti_abAA);
            RAA->to_matrix(RL2AA);
            RAA.reset();
            // TAA
            LAA = std::make_shared<Tensor2d>("New L2 <IJ|AB>", ntri_anti_ijAA, ntri_anti_abAA);
            LAA->read(psio_, PSIF_DFOCC_AMPS);
            auto L2AA = std::make_shared<Matrix>("L2AA", ntri_anti_ijAA, ntri_anti_abAA);
            LAA->to_matrix(L2AA);

            // RBB
            RBB = std::make_shared<Tensor2d>("RL2 <ij|ab>", ntri_anti_ijBB, ntri_anti_abBB);
            RBB->read(psio_, PSIF_DFOCC_AMPS);
            auto RL2BB = std::make_shared<Matrix>("RL2BB", ntri_anti_ijBB, ntri_anti_abBB);
            RBB->to_matrix(RL2BB);
            RBB.reset();
            // TBB
            LBB = std::make_shared<Tensor2d>("New L2 <ij|ab>", ntri_anti_ijBB, ntri_anti_abBB);
            LBB->read(psio_, PSIF_DFOCC_AMPS);
            auto L2BB = std::make_shared<Matrix>("L2BB", ntri_anti_ijBB, ntri_anti_abBB);
            LBB->to_matrix(L2BB);

            // RAB
            RAB = std::make_shared<Tensor2d>("RL2 <Ij|Ab>", naoccA * naoccB, navirA * navirB);
            RAB->read(psio_, PSIF_DFOCC_AMPS);
            auto RL2AB = std::make_shared<Matrix>("RL2AB", naoccA * naoccB, navirA * navirB);
            RAB->to_matrix(RL2AB);
            RAB.reset();
            // TAB
            LAB = std::make_shared<Tensor2d>("New L2 <Ij|Ab>", naoccA * naoccB, navirA * navirB);
            LAB->read(psio_, PSIF_DFOCC_AMPS);
            auto L2AB = std::make_shared<Matrix>("L2AB", naoccA * naoccB, navirA * navirB);
            LAB->to_matrix(L2AB);

            // add entry
            //if (do_diis_ == 1)
                ccsdlDiisManager->add_entry(RL2AA.get(), RL2BB.get(), RL2AB.get(), L2AA.get(), L2BB.get(), L2AB.get());
            RL2AA.reset();
            RL2BB.reset();
            RL2BB.reset();

            // extrapolate
            //if (do_diis_ == 1) {
                if (ccsdlDiisManager->subspace_size() >= cc_mindiis_)
                    ccsdlDiisManager->extrapolate(L2AA.get(), L2BB.get(), L2AB.get());
                LAA->set2(L2AA);
                LBB->set2(L2BB);
                LAB->set2(L2AB);

            //}
            L2AA.reset();
            L2BB.reset();
            L2AB.reset();

            LAA->write(psio_, PSIF_DFOCC_AMPS);
            LBB->write(psio_, PSIF_DFOCC_AMPS);
            LAB->write(psio_, PSIF_DFOCC_AMPS);
            LAA.reset();
            LBB.reset();
            LAB.reset();
    }// end diis

    //=========================
    // Reset & Energy
    //=========================
    ccd_lambda_energy();

} // end ccdl_l2AB_amps

void DFOCC::ccd_lambda_energy() {
    SharedTensor2d K, L, M, I, T, Tnew, T1, T2, U, Tau, W, X, Y;
    SharedTensor2d R, RAA, RBB, RAB, TAA, TBB, TAB;

        // Reset T2AA
        Tnew = std::make_shared<Tensor2d>("New L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        rms_t2AA = Tnew->rms(T);
        T->copy(Tnew);
        Tnew.reset();
        T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

        // AA Energy
        L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
        L->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
        M = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA);
        M->sort(1324, L, 1.0, 0.0);
        L.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ||AB>", naoccA, naoccA, navirA, navirA);
        tei_pqrs_anti_symm_direct(K, M);
        M.reset();
        EccdAA = 0.25 * T->vector_dot(K);
        K.reset();

        // Form L2(IA|JB)
        U = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        U->sort(1324, T, 1.0, 0.0);
        T.reset();
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();

        // Reset T2BB
        Tnew = std::make_shared<Tensor2d>("New L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        Tnew->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
        rms_t2BB = Tnew->rms(T);
        T->copy(Tnew);
        Tnew.reset();
        T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);

        // BB Energy
        L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ia|jb)", naoccB, navirB, naoccB, navirB);
        L->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
        M = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <ij|ab>", naoccB, naoccB, navirB, navirB);
        M->sort(1324, L, 1.0, 0.0);
        L.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <ij||ab>", naoccB, naoccB, navirB, navirB);
        tei_pqrs_anti_symm_direct(K, M);
        M.reset();
        EccdBB = 0.25 * T->vector_dot(K);
        K.reset();

        // Form L2(ia|jb)
        U = std::make_shared<Tensor2d>("L2 (ia|jb)", naoccB, navirB, naoccB, navirB);
        U->sort(1324, T, 1.0, 0.0);
        T.reset();
        U->write_symm(psio_, PSIF_DFOCC_AMPS);
        U.reset();

        // Reset L2AB
        Tnew = std::make_shared<Tensor2d>("New L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        Tnew->read(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_AMPS);
        rms_t2AB = Tnew->rms(T);
        T->copy(Tnew);
        Tnew.reset();
        T->write(psio_, PSIF_DFOCC_AMPS);

        // AB Energy
        L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|jb)", naoccA, navirA, naoccB, navirB);
        L->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
        K = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        K->sort(1324, L, 1.0, 0.0);
        L.reset();
        EccdAB = T->vector_dot(K);
        K.reset();

        // Overall energy
        if (orb_opt_ == "FALSE") {
            EcorrL = EccdAA + EccdBB + EccdAB;
            EccdL = Eref + EcorrL;
        }

        // Form T2(IA|jb)
        U = std::make_shared<Tensor2d>("T2 (IA|jb)", naoccA, navirA, naoccB, navirB);
        U->sort(1324, T, 1.0, 0.0);
        T.reset();
        U->write(psio_, PSIF_DFOCC_AMPS);
        U.reset();

        // combined rms
        double rms_ss = MAX0(rms_t2AA, rms_t2BB);
        rms_l2 = MAX0(rms_ss, rms_t2AB);

}  // end ccd_lambda_energy


}  // namespace dfoccwave
}  // namespace psi
