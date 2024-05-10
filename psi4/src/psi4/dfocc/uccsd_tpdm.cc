/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "dfocc.h"
#include "defines.h"

namespace psi{
namespace dfoccwave {

void DFOCC::uccsd_tpdm()
{
    SharedTensor2d T, T2, L2, Tau, X, Y, Z, V, U, L, G, G2, A, K;

    /////////////////////////////////
    //// OO-Block ///////////////////
    /////////////////////////////////

    //// Alpha BLock ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // G_IJ^Q = 0.5 * P+(IJ) 2*\cal(V)_IJ^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|IJ)", nQ, naoccA, naoccA);
    V = std::make_shared<Tensor2d>("calV (Q|IJ)", nQ, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 1.0);
    V.reset();
    // G_IJ^Q += 0.5 * P+(IJ) V_IJ^Q
    V = std::make_shared<Tensor2d>("V (Q|IJ)", nQ, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 0.5);
    V.reset();
    // G_IJ^Q += 0.5 * P+(IJ) Vt_IJ^Q
    V = std::make_shared<Tensor2d>("Vt (Q|IJ)", nQ, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 0.5);
    V.reset();

    // G_IJ^Q -= 0.5 * P+(IJ) 2*V'_IJ^Q
    V = std::make_shared<Tensor2d>("Vp (Q|IJ)", nQ, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -1.0);
    V.reset();
    // G_IJ^Q -= 0.5 * P+(IJ) Z_IJ^Q
    V = std::make_shared<Tensor2d>("Zeta (Q|IJ)", nQ, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -0.5);
    V.reset();
    // G_IJ^Q -= 0.5 * P+(IJ) G_IJ t_Q
    G->dirprd123(T1c, GijA, -0.5, 1.0);

    // G_IJ^Q += 0.5 * P+(IJ) \sum(M) G_MJ t_IM^Q
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract(false, false, nQ * naoccA, naoccA, naoccA, T, GijA, 0.5, 1.0);
    T.reset();

    // G_IJ^Q -= 0.5 * P+(IJ) \sum(E) l_J^E (t_IE^Q + Tau_IE^Q)
    T = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("Tau (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    G->contract(false, true, nQ * naoccA, naoccA, navirA, T, l1A, -0.5, 1.0);
    T.reset();

    // G_IJ^Q -= 0.5 * P+(IJ) \sum(E) t_J^E (L_IE^Q + 2*V_EI^Q)
    V = std::make_shared<Tensor2d>("V (Q|AI)", nQ, navirA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
    T->swap_3index_col(V);
    V.reset();
    T->scale(2.0);
    U = std::make_shared<Tensor2d>("L2 (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->add(U);
    U.reset();
    G->contract(false, true, nQ * naoccA, naoccA, navirA, T, t1A, -0.5, 1.0);
    T.reset();

    // T3 contribution
    if (wfn_type_ == "DF-CCSD(T)") {
        // G(Q,IM) += P+(IM) 0.5 \sum(JA) B(Q,JA)  * M(IJ,AM)
        // G(Q,IM) += P+(IM) 0.5 \sum(JA) B(Q,JA)  * M(JA,IM)
        T = std::make_shared<Tensor2d>("M <IJ|AM>", naoccA, naoccA, navirA, naoccA);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M <JA|IM>", naoccA, navirA, naoccA, naoccA);
        U->sort(2314, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, 0.5, 1.0);
        U.reset();
        K.reset();

        // G(Q,IM) +=  P+(IM) 0.5 \sum(ja) B(Q,ja) * M(Ij,aM)
        // G(Q,IM) += -P+(IM) 0.5 \sum(ja) B(Q,ja) * M(jI,aM)
        // G(Q,IM) += -P+(IM) 0.5 \sum(ja) B(Q,ja) * M(ja,IM)
        T = std::make_shared<Tensor2d>("M <iJ|aM>", naoccB, naoccA, navirB, naoccA);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M <ia|JM>", naoccB, navirB, naoccA, naoccA);
        U->sort(1324, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, -0.5, 1.0);
        U.reset();
        K.reset();
    }

    // SYMMETRIZE
    G->symmetrize3(G);
    G->scale(2.0);
    G2 = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA);
    G2->set3_act_oo(nfrzc, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if (print_ > 3) G2->print();
    //G2->scale(2.0);
    //G2->print();
    G2.reset();

    //// Beta BLock  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // G_ij^Q = 0.5 * P+(ij) 2*\cal(V)_ij^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|ij)", nQ, naoccB, naoccB);
    V = std::make_shared<Tensor2d>("calV (Q|ij)", nQ, naoccB, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 1.0);
    V.reset();
    // G_ij^Q += 0.5 * P+(ij) V_ij^Q
    V = std::make_shared<Tensor2d>("V (Q|ij)", nQ, naoccB, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 0.5);
    V.reset();
    // G_ij^Q += 0.5 * P+(ij) Vt_ij^Q
    V = std::make_shared<Tensor2d>("Vt (Q|ij)", nQ, naoccB, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, 0.5);
    V.reset();

    // G_ij^Q -= 0.5 * P+(ij) 2*V'_ij^Q
    V = std::make_shared<Tensor2d>("Vp (Q|ij)", nQ, naoccB, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -1.0);
    V.reset();
    // G_ij^Q -= 0.5 * P+(ij) Z_ij^Q
    V = std::make_shared<Tensor2d>("Zeta (Q|ij)", nQ, naoccB, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -0.5);
    V.reset();
    // G_ij^Q -= 0.5 * P+(ij) G_ij t_Q
    G->dirprd123(T1c, GijB, -0.5, 1.0);

    // G_ij^Q += 0.5 * P+(ij) \sum(m) G_mj t_im^Q
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract(false, false, nQ * naoccB, naoccB, naoccB, T, GijB, 0.5, 1.0);
    T.reset();

    // G_ij^Q -= 0.5 * P+(ij) \sum(e) l_j^e (t_ie^Q + Tau_ie^Q)
    T = std::make_shared<Tensor2d>("T1 (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("Tau (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    G->contract(false, true, nQ * naoccB, naoccB, navirB, T, l1B, -0.5, 1.0);
    T.reset();

    // G_ij^Q -= 0.5 * P+(ij) \sum(e) t_j^e (L_ie^Q + 2*V_ei^Q)
    V = std::make_shared<Tensor2d>("V (Q|ai)", nQ, navirB, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|ia)", nQ, naoccB, navirB);
    T->swap_3index_col(V);
    V.reset();
    T->scale(2.0);
    U = std::make_shared<Tensor2d>("L2 (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->add(U);
    U.reset();
    G->contract(false, true, nQ * naoccB, naoccB, navirB, T, t1B, -0.5, 1.0);
    T.reset();

    // T3 contribution
    if (wfn_type_ == "DF-CCSD(T)") {
        // G(Q,im) += P+(im) 0.5 \sum(ja) B(Q,ja)  * M(ij,am)
        // G(Q,im) += P+(im) 0.5 \sum(ja) B(Q,ja)  * M(ja,im)
        T = std::make_shared<Tensor2d>("M <ij|am>", naoccB, naoccB, navirB, naoccB);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M <ja|im>", naoccB, navirB, naoccB, naoccB);
        U->sort(2314, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, 0.5, 1.0);
        U.reset();
        K.reset();

        // G(Q,im) +=  P+(im) 0.5 \sum(JA) M(iJ,Am) * B(Q,JA)
        // G(Q,im) += -P+(im) 0.5 \sum(JA) B(Q,JA)  * M(Ji,Am)
        // G(Q,im) += -P+(im) 0.5 \sum(JA) B(Q,JA)  * M(JA,im)
        T = std::make_shared<Tensor2d>("M <Ij|Am>", naoccA, naoccB, navirA, naoccB);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M <IA|jm>", naoccA, navirA, naoccB, naoccB);
        U->sort(1324, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, -0.5, 1.0);
        U.reset();
        K.reset();
    }

    // symmetrize
    G->symmetrize3(G);
    G->scale(2.0);
    G2 = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|oo)", nQ, noccB, noccB);
    G2->set3_act_oo(nfrzc, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if (print_ > 3) G2->print();
    //G2->scale(2.0);
    //G2->print();

    /* Adding OO for RHF verification
    A = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|OO)", nQ, noccA, noccA);
    A->read(psio_, PSIF_DFOCC_DENS);
    A->add(G2);
    A->scale(1.0);
    A->print();
    A.reset();
    */

    G2.reset();

    /////////////////////////////////
    //// OV-Block ///////////////////
    /////////////////////////////////

    //// Alpha BLock ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // G_IA^Q = 0.5 * Tau_IA^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|IA)", nQ, naoccA, navirA);
    T = std::make_shared<Tensor2d>("Tau (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 0.5);
    T.reset();
    //outfile->Printf("\tStep 1-1 \n");
    //G->scale(2.0);
    //G->print();

    // G_IA^Q += 0.5 * L_IA^Q
    T = std::make_shared<Tensor2d>("L2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 0.5);
    T.reset();
    //outfile->Printf("\tStep 1-2 \n");
    //G->scale(2.0);
    //G->print();

    // G_IA^Q += 0.5 * Z_IA^Q
    T = std::make_shared<Tensor2d>("Zeta (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 0.5);
    T.reset();
    //outfile->Printf("\tStep 1-3 \n");
    //G->scale(2.0);
    //G->print();

    // G_IA^Q += 0.5 * 2*y_IA^Q
    T = std::make_shared<Tensor2d>("Y (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    //T->print();
    G->axpy(T, 1.0);
    T.reset();
    //outfile->Printf("\tStep 1-4 \n");
    //G->scale(2.0);
    //G->print();

    // G_IA^Q -= 0.5 * tEta_IA^Q
    T = std::make_shared<Tensor2d>("Eta2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, -0.5);
    T.reset();

    //outfile->Printf("\tStep 1 \n");
    //G->scale(2.0);
    //G->print();

    // G_IA^Q += 0.5 *  l_I^A t_Q
    G->dirprd123(T1c, l1A, 0.5, 1.0);
    // G_IA^Q += 0.5 * t_I^A (l_Q - G_Q - Gt_Q)
    SharedTensor1d TtA = std::make_shared<Tensor1d>("TEMP", nQ);
    TtA->copy(L1c);
    TtA->axpy(gQ, -1.0);
    TtA->axpy(gQt, -1.0);
    G->dirprd123(TtA, t1A, 0.5, 1.0);
    TtA.reset();

    //outfile->Printf("\tStep 2 \n");
    //G->print();

    // G_IA^Q += 0.5 * \sum(M) t_M^A (G_IM^Q + V_IM^Q - 2*Vp_IM^Q - n_MI^Q)
    T = std::make_shared<Tensor2d>("Temp (Q|IJ)", nQ, naoccA, naoccA);
    U = std::make_shared<Tensor2d>("Eta (Q|IJ)", nQ, naoccA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = std::make_shared<Tensor2d>("G (Q|IJ)", nQ, naoccA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("V (Q|IJ)", nQ, naoccA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("Vp (Q|IJ)", nQ, naoccA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -2.0);
    U.reset();
    G->contract(false, false, nQ * naoccA, navirA, naoccA, T, t1A, 0.5, 1.0);
    T.reset();

    //outfile->Printf("\tStep 3 \n");
    //G->print();

    // G_IA^Q -= 0.5 * \sum(M) Tau_MA^Q Gt_IM
    T = std::make_shared<Tensor2d>("Tau (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract233(false, false, naoccA, navirA, GtijA, T, -0.5, 1.0);

    // G_IA^Q += 0.5 * \sum(E) Tau_IE^Q Gt_EA
    G->contract(false, false, nQ * naoccA, navirA, navirA, T, GtabA, 0.5, 1.0);
    T.reset();

    // G_IA^Q += \sum(E) t_I^E V_EA^Q
    T = std::make_shared<Tensor2d>("V (Q|AB)", nQ, navirA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract233(false, false, naoccA, navirA, t1A, T, 1.0, 1.0);
    T.reset();

    // G_IA^Q += 0.5 * \sum(E) t_IE^Q G_EA
    T = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract(false, false, nQ * naoccA, navirA, navirA, T, GabA, 0.5, 1.0);
    T.reset();

    //outfile->Printf("\tStep 4 \n");
    //G->print();

    // G_IA^Q += 0.5 * \sum(ME) (t_ME^Q - t_EM^Q) * {L2(IM,AE) + 2*V(IE,MA)}
    U = std::make_shared<Tensor2d>("T1 (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();

    V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L->sort(1342, V, 2.0, 1.0);
    V.reset();
    Y = std::make_shared<Tensor2d>("X (ME|IA)", naoccA, navirA, naoccA, navirA);
    Y->sort(2413, L, 1.0, 0.0);
    L.reset();
    G->gemm(false, false, T, Y, 0.5, 1.0);
    Y.reset();
    T.reset();

    // G_IA^Q += 0.5 * \sum(me) (t_me^Q - t_em^Q) * {L2(Im,Ae) - 2*V(Ie,mA)}
    U = std::make_shared<Tensor2d>("T1 (Q|ai)", nQ, navirB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|ia)", nQ, naoccB, navirB);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = std::make_shared<Tensor2d>("T1 (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();

    V = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    L->sort(1342, V, 2.0, 1.0);
    V.reset();
    Y = std::make_shared<Tensor2d>("X (me|IA)", naoccB, navirB, naoccA, navirA);
    Y->sort(2413, L, 1.0, 0.0);
    L.reset();
    G->gemm(false, false, T, Y, 0.5, 1.0);
    Y.reset();
    T.reset();

    //outfile->Printf("\tStep 5 \n");
    //G->print();

    // G_IA^Q += 0.5 * \sum(ME) Tau(IM,AE) (Gt_EM^Q - Gt_ME^Q + l_ME^Q - l_EM^Q)
    U = std::make_shared<Tensor2d>("Gt (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("L1 (Q|AI)", nQ, navirA, naoccA);
    X->read(psio_, PSIF_DFOCC_AMPS);
    U->axpy(X, -1.0);
    X.reset();
    T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("Gt (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("L1 (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);
    T2.reset();
    U = std::make_shared<Tensor2d>("Tau (ME|IA)", naoccA, navirA, naoccA, navirA);
    U->sort(2413, Tau, 1.0, 0.0);
    Tau.reset();
    G->gemm(false, false, T, U, 0.5, 1.0);
    U.reset();
    T.reset();

    // G_IA^Q += 0.5 * \sum(me) Tau(Im,Ae) (Gt_em^Q - Gt_me^Q + l_me^Q - l_em^Q)
    U = std::make_shared<Tensor2d>("Gt (Q|ai)", nQ, navirB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("L1 (Q|ai)", nQ, navirB, naoccB);
    X->read(psio_, PSIF_DFOCC_AMPS);
    U->axpy(X, -1.0);
    X.reset();
    T = std::make_shared<Tensor2d>("Temp (Q|ia)", nQ, naoccB, navirB);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("Gt (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("L1 (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    U = std::make_shared<Tensor2d>("Tau (me|IA)", naoccB, navirB, naoccA, navirA);
    U->sort(2413, Tau, 1.0, 0.0);
    Tau.reset();
    G->gemm(false, false, T, U, 0.5, 1.0);
    U.reset();
    T.reset();

    //outfile->Printf("\tStep 6 \n");
    //G->print();

    // T3 contribution
    if (wfn_type_ == "DF-CCSD(T)") {
        // G(Q,IA) += \sum(JM) B(Q,JM)  * M(IJ,AM)
        // G(Q,IA) += \sum(JM) B(Q,JM)  * M(JM,IA)
        T = std::make_shared<Tensor2d>("M <IJ|AM>", naoccA, naoccA, navirA, naoccA);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M <JM|IA>", naoccA, naoccA, naoccA, navirA);
        U->sort(2413, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, -0.5, 1.0);
        U.reset();
        K.reset();

        // G(Q,IA) += \sum(JB) B(Q,JB)  * M(IJ,AB)
        // G(Q,IA) += \sum(JB) B(Q,JB)  * M(JB,IA)
        T = std::make_shared<Tensor2d>("M <IJ||AB>", naoccA, naoccA, navirA, navirA);
        T->read_anti_symm(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M (JB|IA)", naoccA, navirA, naoccA, navirA);
        U->sort(2413, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, 0.5, 1.0);
        U.reset();
        K.reset();

        // G         += \sum(BD) B(Q,BD) * M[I](A,BD)
        // G[I](Q,A) += \sum(BD) B(Q,BD) * M[I](BD,A)
        U = std::make_shared<Tensor2d>("M[I] <A|BD>", navirA, navirA * navirA);
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
        K->read(psio_, PSIF_DFOCC_INTS, true, true);
        G2 = std::make_shared<Tensor2d>("G[I] (Q|A)", nQ, navirA);
        for (int i = 0; i < naoccA; i++) {
            U->myread(psio_, PSIF_DFOCC_MIABC_AAAA, (size_t)(i * navirA * navirA * navirA) * sizeof(double));
            G2->gemm(false, true, K, U, 0.5, 0.0);
            for (int Q = 0; Q < nQ; Q++) {
                for (int a = 0; a < navirA; a++) {
                    int ia = ia_idxAA->get(i, a);
                    G->add(Q, ia, G2->get(Q, a));
                }
            }
        }
        U.reset();
        K.reset();
        G2.reset();

        // G(Q,IA) += \sum(jm) B(Q,jm) * M(Ij,Am)
        // G(Q,IA) += \sum(jm) B(Q,jm) * M(jm,IA)
        T = std::make_shared<Tensor2d>("M <Ij|Am>", naoccA, naoccB, navirA, naoccB);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M <jm|IA>", naoccB, naoccB, naoccA, navirA);
        U->sort(2413, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, -0.5, 1.0);
        U.reset();
        K.reset();

        // G(Q,IA) += \sum(jb) B(Q,jb) * M(Ij,Ab)
        // G(Q,IA) += \sum(jb) B(Q,jb) * M(jb,IA)
        T = std::make_shared<Tensor2d>("M <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M (jb|IA)", naoccB, navirB, naoccA, navirA);
        U->sort(2413, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, 0.5, 1.0);
        U.reset();
        K.reset();

        // G         += \sum(bd) B(Q,bd) * M[I](A,bd)
        // G[I](Q,A) += \sum(bd) B(Q,bd) * M[I](bd,A)
        U = std::make_shared<Tensor2d>("M[I] <A|bd>", navirA, navirB * navirB);
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB);
        K->read(psio_, PSIF_DFOCC_INTS, true, true);
        G2 = std::make_shared<Tensor2d>("G[I] (Q|A)", nQ, navirA);
        for (int i = 0; i < naoccA; i++) {
            U->myread(psio_, PSIF_DFOCC_MIABC_AABB, (size_t)(i * navirA * navirB * navirB) * sizeof(double));
            G2->gemm(false, true, K, U, 0.5, 0.0);
            for (int Q = 0; Q < nQ; Q++) {
                for (int a = 0; a < navirA; a++) {
                    int ia = ia_idxAA->get(i, a);
                    G->add(Q, ia, G2->get(Q, a));
                }
            }
        }
        U.reset();
        K.reset();
        G2.reset();
    }

    // Form overall OV Block
    G2 = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA);
    G2->set3_act_ov(nfrzc, naoccA, navirA, nvirA, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if (print_ > 3) G2->print();
    //G2->scale(2.0);
    //G2->print();

    // Form G_vo^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA);
    G->swap_3index_col(G2);
    G2.reset();
    G->write(psio_, PSIF_DFOCC_DENS);
    if (print_ > 3) G->print();
    G.reset();

    //// Beta BLock  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // G_ia^Q = 0.5 * Tau_ia^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|ia)", nQ, naoccB, navirB);
    T = std::make_shared<Tensor2d>("Tau (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 0.5);
    T.reset();
    // G_ia^Q += 0.5 * L_ia^Q
    T = std::make_shared<Tensor2d>("L2 (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 0.5);
    T.reset();
    // G_ia^Q += 0.5 * Z_ia^Q
    T = std::make_shared<Tensor2d>("Zeta (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 0.5);
    T.reset();
    // G_ia^Q += 0.5 * 2*y_ia^Q
    T = std::make_shared<Tensor2d>("Y (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, 1.0);
    T.reset();
    // G_ia^Q -= 0.5 * tEta_ia^Q
    T = std::make_shared<Tensor2d>("Eta2 (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(T, -0.5);
    T.reset();

    // G_ia^Q += 0.5 *  l_i^a t_Q
    G->dirprd123(T1c, l1B, 0.5, 1.0);
    // G_ia^Q += 0.5 * t_i^a (l_Q - G_Q - Gt_Q)
    SharedTensor1d TtB = std::make_shared<Tensor1d>("TEMP", nQ);
    TtB->copy(L1c);
    TtB->axpy(gQ, -1.0);
    TtB->axpy(gQt, -1.0);
    G->dirprd123(TtB, t1B, 0.5, 1.0);
    TtB.reset();

    // G_ia^Q += 0.5 * \sum(m) t_m^a (G_im^Q + V_im^Q - 2*Vp_im^Q - n_mi^Q)
    T = std::make_shared<Tensor2d>("Temp (Q|ij)", nQ, naoccB, naoccB);
    U = std::make_shared<Tensor2d>("Eta (Q|ij)", nQ, naoccB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = std::make_shared<Tensor2d>("G (Q|ij)", nQ, naoccB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("V (Q|ij)", nQ, naoccB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("Vp (Q|ij)", nQ, naoccB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -2.0);
    U.reset();
    G->contract(false, false, nQ * naoccB, navirB, naoccB, T, t1B, 0.5, 1.0);
    T.reset();

    // G_ia^Q -= 0.5 * \sum(m) Tau_ma^Q Gt_im
    T = std::make_shared<Tensor2d>("Tau (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract233(false, false, naoccB, navirB, GtijB, T, -0.5, 1.0);

    // G_ia^Q += 0.5 * \sum(e) Tau_ie^Q Gt_ea
    G->contract(false, false, nQ * naoccB, navirB, navirB, T, GtabB, 0.5, 1.0);
    T.reset();

    // G_ia^Q += \sum(e) t_i^e V_ea^Q
    T = std::make_shared<Tensor2d>("V (Q|ab)", nQ, navirB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract233(false, false, naoccB, navirB, t1B, T, 1.0, 1.0);
    T.reset();

    // G_ia^Q += 0.5 * \sum(e) t_ie^Q G_ea
    T = std::make_shared<Tensor2d>("T1 (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    G->contract(false, false, nQ * naoccB, navirB, navirB, T, GabB, 0.5, 1.0);
    T.reset();

    // G_ia^Q += 0.5 * \sum(me) (t_me^Q - t_em^Q) * {L2(im,ae) + 2*V(ie,ma)}
    T = std::make_shared<Tensor2d>("Temp (Q|ia)", nQ, naoccB, navirB);
    U = std::make_shared<Tensor2d>("T1 (Q|ai)", nQ, navirB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = std::make_shared<Tensor2d>("T1 (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();

    V = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    L->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    L->sort(1342, V, 2.0, 1.0);
    V.reset();
    Y = std::make_shared<Tensor2d>("X (me|ia)", naoccB, navirB, naoccB, navirB);
    Y->sort(2413, L, 1.0, 0.0);
    L.reset();
    G->gemm(false, false, T, Y, 0.5, 1.0);
    Y.reset();
    T.reset();

    // G_ia^Q += 0.5 * \sum(ME) (t_ME^Q - t_EM^Q) * {L2(Mi,Ea) - 2*V(iE,Ma)}
    T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
    U = std::make_shared<Tensor2d>("T1 (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = std::make_shared<Tensor2d>("T1 (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();

    V = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    L->sort(3124, V, 2.0, 1.0);
    V.reset();
    Y = std::make_shared<Tensor2d>("X (ME|ia)", naoccA, navirA, naoccB, navirB);
    Y->sort(1324, L, 1.0, 0.0);
    L.reset();
    G->gemm(false, false, T, Y, 0.5, 1.0);
    Y.reset();
    T.reset();

    // G_ia^Q += 0.5 * \sum(me) Tau(im,ae) (Gt_em^Q - Gt_me^Q + l_me^Q - l_em^Q)
    U = std::make_shared<Tensor2d>("Gt (Q|ai)", nQ, navirB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("L1 (Q|ai)", nQ, navirB, naoccB);
    X->read(psio_, PSIF_DFOCC_AMPS);
    U->axpy(X, -1.0);
    X.reset();
    T = std::make_shared<Tensor2d>("Temp (Q|ia)", nQ, naoccB, navirB);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("Gt (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("L1 (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
    T2.reset();
    U = std::make_shared<Tensor2d>("Tau (me|ia)", naoccB, navirB, naoccB, navirB);
    U->sort(2413, Tau, 1.0, 0.0);
    Tau.reset();
    G->gemm(false, false, T, U, 0.5, 1.0);
    U.reset();
    T.reset();

    // G_ia^Q += 0.5 * \sum(ME) Tau(Mi,Ea) (Gt_EM^Q - Gt_ME^Q + l_ME^Q - l_EM^Q)
    U = std::make_shared<Tensor2d>("Gt (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("L1 (Q|AI)", nQ, navirA, naoccA);
    X->read(psio_, PSIF_DFOCC_AMPS);
    U->axpy(X, -1.0);
    X.reset();
    T = std::make_shared<Tensor2d>("Temp (Q|IA)", nQ, naoccA, navirA);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("Gt (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, -1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("L1 (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    U = std::make_shared<Tensor2d>("Tau (ME|ia)", naoccA, navirA, naoccB, navirB);
    U->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    G->gemm(false, false, T, U, 0.5, 1.0);
    U.reset();
    T.reset();

    // T3 contribution
    if (wfn_type_ == "DF-CCSD(T)") {
        // G(Q,ia) += \sum(JM) B(Q,JM)  * M(iJ,aM)
        // G(Q,ia) += \sum(JM) B(Q,JM)  * M(JM,ia)
        T = std::make_shared<Tensor2d>("M <iJ|aM>", naoccB, naoccA, navirB, naoccA);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M <JM|ia>", naoccA, naoccA, naoccB, navirB);
        U->sort(2413, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, -0.5, 1.0);
        U.reset();
        K.reset();

        // G(Q,ia) += \sum(JB) B(Q,JB)  * M(iJ,aB)
        // G(Q,ia) += \sum(JB) B(Q,JB)  * M(Ji,Ba)
        // G(Q,ia) += \sum(JB) B(Q,JB)  * M(JB,ia)
        T = std::make_shared<Tensor2d>("M <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M (IA|jb)", naoccA, navirA, naoccB, navirB);
        U->sort(1324, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, 0.5, 1.0);
        U.reset();
        K.reset();

        // G[i](Q,a) += \sum(BD) B(Q,BD)    * M[i](BD,a)
        U = std::make_shared<Tensor2d>("M[I] <a|BD>", navirB, navirA * navirA);
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
        K->read(psio_, PSIF_DFOCC_INTS, true, true);
        G2 = std::make_shared<Tensor2d>("G[I] (Q|a)", nQ, navirB);
        for (int i = 0; i < naoccB; i++) {
            U->myread(psio_, PSIF_DFOCC_MIABC_BBAA, (size_t)(i * navirB * navirA * navirA) * sizeof(double));
            G2->gemm(false, true, K, U, 0.5, 0.0);
            for (int Q = 0; Q < nQ; Q++) {
                for (int a = 0; a < navirB; a++) {
                    int ia = ia_idxBB->get(i, a);
                    G->add(Q, ia, G2->get(Q, a));
                }
            }
        }
        U.reset();
        K.reset();
        G2.reset();

        // G(Q,ia) += \sum(jm) B(Q,jm)  * M(ij,am)
        // G(Q,ia) += \sum(jm) B(Q,jm)  * M(jm,ia)
        T = std::make_shared<Tensor2d>("M <ij|am>", naoccB, naoccB, navirB, naoccB);
        T->read(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M <jm|ia>", naoccB, naoccB, naoccB, navirB);
        U->sort(2413, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, -0.5, 1.0);
        U.reset();
        K.reset();

        // G(Q,ia) += \sum(jb) B(Q,jb)  * M(ij,ab)
        // G(Q,ia) += \sum(jb) B(Q,jb)  * M(jb,ia)
        T = std::make_shared<Tensor2d>("M <ij||ab>", naoccB, naoccB, navirB, navirB);
        T->read_anti_symm(psio_, PSIF_DFOCC_DENS);
        U = std::make_shared<Tensor2d>("M (jb|ia)", naoccB, navirB, naoccB, navirB);
        U->sort(2413, T, 1.0, 0.0);
        T.reset();
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB);
        K->read(psio_, PSIF_DFOCC_INTS);
        G->gemm(false, false, K, U, 0.5, 1.0);
        U.reset();
        K.reset();

        // G         += \sum(bd) B(Q,bd) * M[i](a,bd)
        // G[i](Q,a) += \sum(bd) B(Q,bd) * M[i](bd,a)
        U = std::make_shared<Tensor2d>("M[i] <a|bd>", navirB, navirB * navirB);
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB);
        K->read(psio_, PSIF_DFOCC_INTS, true, true);
        G2 = std::make_shared<Tensor2d>("G[i] (Q|a)", nQ, navirB);
        for (int i = 0; i < naoccB; i++) {
            U->myread(psio_, PSIF_DFOCC_MIABC_BBBB, (size_t)(i * navirB * navirB * navirB) * sizeof(double));
            G2->gemm(false, true, K, U, 0.5, 0.0);
            for (int Q = 0; Q < nQ; Q++) {
                for (int a = 0; a < navirB; a++) {
                    int ia = ia_idxBB->get(i, a);
                    G->add(Q, ia, G2->get(Q, a));
                }
            }
        }
        U.reset();
        K.reset();
        G2.reset();
    }

    // Form overall OV Block
    G2 = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|ov)", nQ, noccB, nvirB);
    G2->set3_act_ov(nfrzc, naoccB, navirB, nvirB, G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS);
    if (print_ > 3) G2->print();

    // Form G_vo^Q
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|vo)", nQ, nvirB, noccB);
    G->swap_3index_col(G2);

    /* Adding OV for RHF verification
    SharedTensor2d G2A = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|OV)", nQ, noccA, nvirA);
    G2A->read(psio_, PSIF_DFOCC_DENS);
    G2A->add(G2);
    G2A->print();
    G2A.reset();
    */
    //G2->scale(2.0);
    //G2->print();
    G2.reset();

    /* Adding VO for RHF verification
    SharedTensor2d GA = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|VO)", nQ, nvirA, noccA);
    GA->read(psio_, PSIF_DFOCC_DENS);
    GA->add(G);
    GA->print();
    GA.reset();
    */
    //G->scale(2.0);
    //G->print();
    G->write(psio_, PSIF_DFOCC_DENS);
    if (print_ > 3) G->print();
    G.reset();

    /////////////////////////////////
    //// VV-Block ///////////////////
    /////////////////////////////////

    //// Alpha BLock ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //G{AB}^{Q} =  - 0.5 *  P_{+}(AB) [2 V_{AB}^{Q} + 2 {Vt}_{AB}^{Q} + etatilde_{AB}^{Q} + \cal G_{AB} t_{Q}]
    // G_ab^Q -= 0.5 * P+(ab) 2*V_ab^Q
    V = std::make_shared<Tensor2d>("V (Q|AB)", nQ, navirA, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|AB)", nQ, navirA, navirA);
    G->axpy(V, -1.0);
    V.reset();
    // G_ab^Q -= 0.5 * P+(ab) 2*Vt_ab^Q
    V = std::make_shared<Tensor2d>("Vt (Q|AB)", nQ, navirA, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -1.0);
    V.reset();
    // G_ab^Q -= 0.5 * P+(ab) tEta_ab^Q
    V = std::make_shared<Tensor2d>("Eta2 (Q|AB)", nQ, navirA, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -0.5);
    V.reset();
    // G_ab^Q -= 0.5 * P+(ab) G_ab t_Q
    G->dirprd123(T1c, GabA, -0.5, 1.0);

    //G{AB}^{Q} =  + 0.5 * P_{+}(AB) \sum_{M}^{occ} l_{B}^{M} [tau_{MA}^{Q} - t_{AM}^{Q}]
    U = std::make_shared<Tensor2d>("Tau (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|AI)", nQ, navirA, naoccA);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = std::make_shared<Tensor2d>("T1 (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    G->contract(false, false, nQ * navirA, navirA, naoccA, T, l1A, -0.5, 1.0);
    T.reset();

    //G{AB}^{Q} =  + 0.5 * P_{+}(AB) \sum_{M}^{occ} t_{M}^{B} [eta_{MA}^{Q} + {\cal G}_{AM}^{Q} + L_{MA}^{Q} + 2 V_{AM}^{Q}]
    U = std::make_shared<Tensor2d>("Eta (Q|IA)", nQ, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("L2 (Q|IA)", nQ, naoccA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    U->add(T);
    T.reset();
    T = std::make_shared<Tensor2d>("Temp (Q|AI)", nQ, navirA, naoccA);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("G (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("V (Q|AI)", nQ, navirA, naoccA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 2.0);
    U.reset();
    G->contract(false, false, nQ * navirA, navirA, naoccA, T, t1A, 0.5, 1.0);
    T.reset();

    // PPL
    ccd_tpdm_pplA(G, "Tau");

    // T3 contribution
    if (wfn_type_ == "DF-CCSD(T)") {
        // G(Q,AB) += P+(AB) 0.5 \sum(IC) B(Q,IC) * M(IC,AB)
        U = std::make_shared<Tensor2d>("M[I] <A|BD>", navirA, navirA * navirA);
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_INTS);
        T = std::make_shared<Tensor2d>("B[I] (Q|A)", nQ, navirA);
        for (int i = 0; i < naoccA; i++) {
            U->myread(psio_, PSIF_DFOCC_MIABC_AAAA, (size_t)(i * navirA * navirA * navirA) * sizeof(double));
            for (int Q = 0; Q < nQ; Q++) {
                for (int c = 0; c < navirA; c++) {
                    int ic = ia_idxAA->get(i, c);
                    T->set(Q, c, K->get(Q, ic));
                }
            }
            G->gemm(false, false, T, U, 0.5, 1.0);
        }
        U.reset();
        K.reset();
        T.reset();

        // G(Q,AB) += P+(AB) 0.5 \sum(ic) B(Q,ic) * M(ic,AB)
        U = std::make_shared<Tensor2d>("M[i] <a|BD>", navirB, navirA * navirA);
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB);
        K->read(psio_, PSIF_DFOCC_INTS);
        T = std::make_shared<Tensor2d>("B[i] (Q|a)", nQ, navirB);
        for (int i = 0; i < naoccB; i++) {
            U->myread(psio_, PSIF_DFOCC_MIABC_BBAA, (size_t)(i * navirB * navirA * navirA) * sizeof(double));
            for (int Q = 0; Q < nQ; Q++) {
                for (int c = 0; c < navirB; c++) {
                    int ic = ia_idxBB->get(i, c);
                    T->set(Q, c, K->get(Q, ic));
                }
            }
            G->gemm(false, false, T, U, 0.5, 1.0);
        }
        U.reset();
        K.reset();
        T.reset();
    }

    // SYMMETRIZE
    G->symmetrize3(G);
    G->scale(2.0);
    G2 = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA);
    G2->set3_act_vv(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS, true, true);
    if (print_ > 3) G2->print();
    //G2->scale(2.0);
    //G2->print();
    G2.reset();

    //// Beta BLock  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //G{ab}^{Q} = - 0.5 * P_{+}(ab) [2 V_{ab}^{Q} + 2 {Vt}_{ab}^{Q} + etatilde_{ab}^{Q} + \cal G_{ab} t_{Q}]
    // G_ab^Q -= 0.5 * P+(ab) 2*V_ab^Q
    V = std::make_shared<Tensor2d>("V (Q|ab)", nQ, navirB, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|ab)", nQ, navirB, navirB);
    G->axpy(V, -1.0);
    V.reset();
    // G_ab^Q -= 0.5 * P+(ab) 2*Vt_ab^Q
    V = std::make_shared<Tensor2d>("Vt (Q|ab)", nQ, navirB, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -1.0);
    V.reset();
    // G_ab^Q -= 0.5 * P+(ab) tEta_ab^Q
    V = std::make_shared<Tensor2d>("Eta2 (Q|ab)", nQ, navirB, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    G->axpy(V, -0.5);
    V.reset();
    // G_ab^Q -= 0.5 * P+(ab) G_ab t_Q
    G->dirprd123(T1c, GabB, -0.5, 1.0);

    //G{ab}^{Q} = + 0.5 * P_{+}(ab) \sum_{m}^{occ} l_{b}^{m} [tau_{ma}^{Q} - t_{am}^{Q}]
    U = std::make_shared<Tensor2d>("Tau (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("Temp (Q|ai)", nQ, navirB, naoccB);
    T->swap_3index_col(U);
    U.reset();
    T->scale(-1.0);
    U = std::make_shared<Tensor2d>("T1 (Q|ai)", nQ, navirB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    G->contract(false, false, nQ * navirB, navirB, naoccB, T, l1B, -0.5, 1.0);
    T.reset();

    //G{ab}^{Q} = + 0.5 * P_{+}(ab) \sum_{m}^{occ} t_{m}^{b} [eta_{ma}^{Q} + \cal G_{am}^{Q} + L_{ma}^{Q} + 2 V_{am}^{Q}]
    U = std::make_shared<Tensor2d>("Eta (Q|ia)", nQ, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("L2 (Q|ia)", nQ, naoccB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    U->add(T);
    T.reset();
    T = std::make_shared<Tensor2d>("Temp (Q|ai)", nQ, navirB, naoccB);
    T->swap_3index_col(U);
    U.reset();
    U = std::make_shared<Tensor2d>("G (Q|ai)", nQ, navirB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 1.0);
    U.reset();
    U = std::make_shared<Tensor2d>("V (Q|ai)", nQ, navirB, naoccB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    T->axpy(U, 2.0);
    U.reset();
    G->contract(false, false, nQ * navirB, navirB, naoccB, T, t1B, 0.5, 1.0);
    T.reset();

    // PPL
    ccd_tpdm_pplB(G, "Tau");

    // T3 contribution
    if (wfn_type_ == "DF-CCSD(T)") {
        // G(Q,ab) += P+(ab) 0.5 \sum(ic) B(Q,ic) * M(ic,ab)
        U = std::make_shared<Tensor2d>("M[i] <a|bd>", navirB, navirB * navirB);
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB);
        K->read(psio_, PSIF_DFOCC_INTS);
        T = std::make_shared<Tensor2d>("B[i] (Q|a)", nQ, navirB);
        for (int i = 0; i < naoccB; i++) {
            U->myread(psio_, PSIF_DFOCC_MIABC_BBBB, (size_t)(i * navirB * navirB * navirB) * sizeof(double));
            for (int Q = 0; Q < nQ; Q++) {
                for (int c = 0; c < navirB; c++) {
                    int ic = ia_idxBB->get(i, c);
                    T->set(Q, c, K->get(Q, ic));
                }
            }
            G->gemm(false, false, T, U, 0.5, 1.0);
        }
        U.reset();
        K.reset();
        T.reset();

        // G(Q,ab) += P+(ab) 0.5 \sum(IC) B(Q,IC) * M(IC,ab)
        U = std::make_shared<Tensor2d>("M[I] <A|bd>", navirA, navirB * navirB);
        K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
        K->read(psio_, PSIF_DFOCC_INTS);
        T = std::make_shared<Tensor2d>("B[I] (Q|A)", nQ, navirA);
        for (int i = 0; i < naoccA; i++) {
            U->myread(psio_, PSIF_DFOCC_MIABC_AABB, (size_t)(i * navirA * navirB * navirB) * sizeof(double));
            for (int Q = 0; Q < nQ; Q++) {
                for (int c = 0; c < navirA; c++) {
                    int ic = ia_idxAA->get(i, c);
                    T->set(Q, c, K->get(Q, ic));
                }
            }
            G->gemm(false, false, T, U, 0.5, 1.0);
        }
        U.reset();
        K.reset();
        T.reset();
    }

    // symmetrize
    G->symmetrize3(G);
    G->scale(2.0);
    G2 = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|vv)", nQ, nvirB, nvirB);
    G2->set3_act_vv(G);
    G.reset();
    G2->write(psio_, PSIF_DFOCC_DENS, true, true);
    if (print_ > 3) G2->print();
    //G2->scale(2.0);
    //G2->print();
    /* Adding VV for rhf verification
    SharedTensor2d GA = std::make_shared<Tensor2d>("Correlation 3-Index TPDM (Q|VV)", nQ, nvirA, nvirA);
    GA->read(psio_, PSIF_DFOCC_DENS, true, true);
    GA->add(G2);
    GA->print();
    GA.reset();
    */

    G2.reset();

    // remove files
    remove_binary_file(PSIF_DFOCC_MIABC_AAAA);
    remove_binary_file(PSIF_DFOCC_MIABC_BBBB);
    remove_binary_file(PSIF_DFOCC_MIABC_AABB);
    remove_binary_file(PSIF_DFOCC_MIABC_BBAA);

    //outfile->Printf("\tI am here.\n");
} // end uccsd_tpdm
}  // namespace dfoccwave
}  // namespace psi
