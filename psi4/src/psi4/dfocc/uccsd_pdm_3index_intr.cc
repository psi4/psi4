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

namespace psi{
namespace dfoccwave {

void DFOCC::uccsd_pdm_3index_intr()
{

    SharedTensor2d T, T2, L2, Tau, X, Y, Z, V, Vt, U, L, Lt, Vab, Vtab, Vij;
    // Build Gt_mi    (C1)
    // Gt_MI = G_MI + \sum_{E} t_M^E l_I^E         (25)
    GtijA->copy(GijA);
    GtijA->gemm(false, true, t1A, l1A, 1.0, 1.0);
    //GtijA->print();
    // Gt_mi = G_mi + \sum_{e} t_m^e l_i^e         (26)
    GtijB->copy(GijB);
    GtijB->gemm(false, true, t1B, l1B, 1.0, 1.0);
    //GtijB->print();
    // Build Gt_ae     (C2)
    // Gt_AE = G_AE - \sum_{M} t_M^E l_M^A         (28)
    GtabA->copy(GabA);
    GtabA->gemm(true, false, l1A, t1A, -1.0, 1.0);
    //GtabA->print();
    // Gt_ae = G_ae - \sum_{m} t_m^e l_m^a         (29)
    GtabB->copy(GabB);
    GtabB->gemm(true, false, l1B, t1B, -1.0, 1.0);
    //GtabB->print();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // V(AB,CD) Intermediates /////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    // V(AB,CD) = 1/2 \sum_{M,N} Tau (MN,CD) * L(MN,AB)       (31)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);
    T2.reset();
    L2 = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("V <AB|CD>", navirA, navirA, navirA, navirA);
    V->gemm(true, false, L2, Tau, 0.5, 0.0);
    L2.reset();
    Tau.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset();
    // V(ab,cd) = 1/2 \sum_{m,n} Tau (mn,cd) * L(mn,ab)       (32) 
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
    T2.reset();
    L2 = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    L2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("V <ab|cd>", navirB, navirB, navirB, navirB);
    V->gemm(true, false, L2, Tau, 0.5, 0.0);
    L2.reset();
    Tau.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset();
    // V(Ab,Cd) = \sum_{M,n} Tau (Mn,Cd) * L(Mn,Ab)        (33)
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L2->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("V <Ab|Cd>", navirA, navirB, navirA, navirB);
    V->gemm(true, false, L2, Tau, 0.5, 0.0);
    L2.reset();   // used in V(aB,cD) tensor
    Tau.reset();  // used in V(aB,cD) tensor
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset();
    // V(aB,cD) = \sum_{m,N} Tau (Nm,Dc) * L(Nm,Ba)        (34)
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    L2 = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    L2->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <Ba|Dc>", navirA, navirB, navirA, navirB);
    X->gemm(true, false, L2, Tau, 1.0, 0.0);
    L2.reset();
    Tau.reset();
    V = std::make_shared<Tensor2d>("V <aB|cD>", navirB, navirA, navirB, navirA);
    V->sort(2143, X, 1.0, 0.0);
    X.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset();
*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Vtilde(IA,JB) Intermediates ///////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // V_tilde(IA,JB) = V(IA,JB) + [t(I,B) * l(J,A)]       (36)
    V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Vt = std::make_shared<Tensor2d>("Vt (IA|JB)", naoccA, navirA, naoccA, navirA);
    Vt->copy(V);
    V.reset(); 
    X = std::make_shared<Tensor2d>("X (IB|JA)", naoccA, navirA, naoccA, navirA);
    X->dirprd224(t1A, l1A, 1.0, 0.0);
    Vt->sort(1432, X, 1.0, 1.0);
    X.reset();
    Vt->write(psio_, PSIF_DFOCC_AMPS);
    Vt.reset();
    // V_tilde(ia,jb) = V(ia,jb) + [t(i,b) * t(l,a)]       (37)
    V = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Vt = std::make_shared<Tensor2d>("Vt (ia|jb)", naoccB, navirB, naoccB, navirB);
    Vt->copy(V);
    V.reset();
    X = std::make_shared<Tensor2d>("X (ib|ja)", naoccB, navirB, naoccB, navirB);
    X->dirprd224(t1B, l1B, 1.0, 0.0);
    Vt->sort(1432, X, 1.0, 1.0);
    X.reset();
    Vt->write(psio_, PSIF_DFOCC_AMPS);
    Vt.reset();
    // V_tilde(Ia,Jb) = V(Ia,Jb)                            (38)
    V = std::make_shared<Tensor2d>("V (Ia|Jb)", naoccA, navirB, naoccA, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Vt = std::make_shared<Tensor2d>("Vt (Ia|Jb)", naoccA, navirB, naoccA, navirB);
    Vt->copy(V);
    V.reset();
    Vt->write(psio_, PSIF_DFOCC_AMPS);
    Vt.reset();
    // V_tilde(iA,jB) = V(iA,jB)                            (39)
    V = std::make_shared<Tensor2d>("V (iA|jB)", naoccB, navirA, naoccB, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Vt = std::make_shared<Tensor2d>("Vt (iA|jB)", naoccB, navirA, naoccB, navirA);
    Vt->copy(V);
    V.reset();
    Vt->write(psio_, PSIF_DFOCC_AMPS);
    Vt.reset();
    // V_tilde(Ia,jB) = V(Ia,jB) + [t(I,B) * l(j,a)]         (40)
    V = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Vt = std::make_shared<Tensor2d>("Vt (Ia|jB)", naoccA, navirB, naoccB, navirA);
    Vt->copy(V);
    V.reset();
    X = std::make_shared<Tensor2d>("X (IB|ja)", naoccA, navirA, naoccB, navirB);
    X->dirprd224(t1A, l1B, 1.0, 0.0);
    Vt->sort(1432, X, 1.0, 1.0);
    X.reset();
    Vt->write(psio_, PSIF_DFOCC_AMPS);
    Vt.reset();
    // V_tilde(iA,Jb) = V(iA,Jb) + [t(i,b) * l(J,A)]          (41)
    V = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Vt = std::make_shared<Tensor2d>("Vt (iA|Jb)", naoccB, navirA, naoccA, navirB);
    Vt->copy(V);
    V.reset();
    X = std::make_shared<Tensor2d>("X (ib|JA)", naoccB, navirB, naoccA, navirA);
    X->dirprd224(t1B, l1A, 1.0, 0.0);
    Vt->sort(1432, X, 1.0, 1.0);
    X.reset();
    Vt->write(psio_, PSIF_DFOCC_AMPS);
    Vt.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // V(IJ,KA) Intermediates ////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // V(IJ,KA) = \sum_{E} t(J,E) * V(IE,KA)          (43)
    U = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (EI|KA)", navirA, naoccA, naoccA, navirA);
    X->sort(2134, U, 1.0, 0.0);
    U.reset();
    Y = std::make_shared<Tensor2d>("Y <JI|KA>", naoccA, naoccA, naoccA, navirA);
    Y->contract(false, false, naoccA, naoccA * naoccA * navirA, navirA, t1A, X, 1.0, 0.0);
    X.reset(); 
    V = std::make_shared<Tensor2d>("V <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    V->sort(2134, Y, 1.0, 0.0);
    Y.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset();
    // V(ij,ka) = \sum_{e} t(j,e) * V(ie,ka)          (44)
    U = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ei|ka)", navirB, naoccB, naoccB, navirB);
    X->sort(2134, U, 1.0, 0.0);
    U.reset();
    Y = std::make_shared<Tensor2d>("Y <ji|ka>", naoccB, naoccB, naoccB, navirB);
    Y->contract(false, false, naoccB, naoccB * naoccB * navirB, navirB, t1B, X, 1.0, 0.0);
    X.reset();
    V = std::make_shared<Tensor2d>("V <ij|ka>", naoccB, naoccB, naoccB, navirB);
    V->sort(2134, Y, 1.0, 0.0);
    Y.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset();
    // V(Ij,Ka) = \sum_{E} t(j,e) * V(Ie,Ka)         (45)
    U = std::make_shared<Tensor2d>("V (Ia|Jb)", naoccA, navirB, naoccA, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (eI|Ka)", navirB, naoccA, naoccA, navirB);
    X->sort(2134, U, 1.0, 0.0);
    U.reset();
    Y = std::make_shared<Tensor2d>("Y <jI|Ka>", naoccB, naoccA, naoccA, navirB);
    Y->contract(false, false, naoccB, naoccA * naoccA * navirB, navirB, t1B, X, 1.0, 0.0);
    X.reset();
    V = std::make_shared<Tensor2d>("V <Ij|Ka>", naoccA, naoccB, naoccA, navirB);
    V->sort(2134, Y, 1.0, 0.0);
    Y.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset();
    // V(iJ,kA) = \sum_{E} t(J,E) * V(iE,kA)         (46)
    U = std::make_shared<Tensor2d>("V (iA|jB)", naoccB, navirA, naoccB, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (Ei|kA)", navirA, naoccB, naoccB, navirA);
    X->sort(2134, U, 1.0, 0.0);
    U.reset();
    Y = std::make_shared<Tensor2d>("Y (Ji|kA)", naoccA, naoccB, naoccB, navirA);
    Y->contract(false, false, naoccA, naoccB * naoccB * navirA, navirA, t1A, X, 1.0, 0.0);
    X.reset();
    V = std::make_shared<Tensor2d>("V <iJ|kA>", naoccB, naoccA, naoccB, navirA);
    V->sort(2134, Y, 1.0, 0.0);
    Y.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset();
    // V(Ij,kA) = \sum_{e} t(j,e) * V(Ie,kA)         (47)
    U = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (eI|kA)", navirB, naoccA, naoccB, navirA);
    X->sort(2134, U, 1.0, 0.0);
    U.reset();
    Y = std::make_shared<Tensor2d>("Y (jI|kA)", naoccB, naoccA, naoccB, navirA);
    Y->contract(false, false, naoccB, naoccA * naoccB * navirA, navirB, t1B, X, 1.0, 0.0);
    X.reset();
    V = std::make_shared<Tensor2d>("V <Ij|kA>", naoccA, naoccB, naoccB, navirA);
    V->sort(2134, Y, 1.0, 0.0);
    Y.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset(); 
    // V(iJ,Ka) = \sum_{E} t(J,E) * V(iE,Ka)          (48)
    U = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (Ei|Ka)", navirA, naoccB, naoccA, navirB);
    X->sort(2134, U, 1.0, 0.0);
    U.reset(); 
    Y = std::make_shared<Tensor2d>("Y (Ji|Ka)", naoccA, naoccB, naoccA, navirB);
    Y->contract(false, false, naoccA, naoccB * naoccA * navirB, navirA, t1A, X, 1.0, 0.0);
    X.reset();
    V = std::make_shared<Tensor2d>("V <iJ|Ka>", naoccB, naoccA, naoccA, navirB);
    V->sort(2134, Y, 1.0, 0.0);
    Y.reset();
    V->write(psio_, PSIF_DFOCC_AMPS);
    V.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Lambda_tilde (IE,AF)  Intermediates ///////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Lt_IEAF = 0.5 * \sum{M,N} Tau(MN,AF) * L(MN,IE) (50)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);
    T2.reset();
    L = std::make_shared<Tensor2d>("L <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    L->read(psio_, PSIF_DFOCC_AMPS);
    Lt = std::make_shared<Tensor2d>("Lambda <IE|AF>", naoccA, navirA, navirA, navirA);
    Lt->gemm(true, false, L, Tau, 0.5, 0.0);
    L.reset();
    Tau.reset();
    Lt->write(psio_, PSIF_DFOCC_AMPS);
    Lt.reset();

    // Lt_ieaf = 0.5 * \sum{m,n} Tau(mn,af) * L(mn,ie) (51)
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
    T2.reset();
    L = std::make_shared<Tensor2d>("L <ij|ka>", naoccB, naoccB, naoccB, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    Lt = std::make_shared<Tensor2d>("Lambda <ie|af>", naoccB, navirB, navirB, navirB);
    Lt->gemm(true, false, L, Tau, 0.5, 0.0);
    L.reset();
    Tau.reset();
    Lt->write(psio_, PSIF_DFOCC_AMPS);
    //Lt->print();
    Lt.reset();

    // Lt_IeAf = \sum{M,n} Tau(Mn,AF) * L(Mn,Ie)       (52)
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    L = std::make_shared<Tensor2d>("L <Ij|Ka>", naoccA, naoccB, naoccA, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    Lt = std::make_shared<Tensor2d>("Lambda <Ie|Af>", naoccA, navirB, navirA, navirB);
    Lt->gemm(true, false, L, Tau, 1.0, 0.0);
    L.reset();
    //Tau.reset(); // Asagida kullaniyorum
    Lt->write(psio_, PSIF_DFOCC_AMPS);
    //Lt->print();
    Lt.reset();

    // Lt_iEaF = \sum{M,n} Tau(Mn,Fa) * L(nM,iE)       (53)
    X = std::make_shared<Tensor2d>("L <iJ|kA>", naoccB, naoccA, naoccB, navirA);
    X->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("X <iE|Mn>", naoccB, navirA, naoccA, naoccB);
    L->sort(3421, X, 1.0, 0.0);
    X.reset();
    Y = std::make_shared<Tensor2d>("Y <iE|Fa>", naoccB, navirA, navirA, navirB);
    Y->gemm(false, false, L, Tau, 1.0, 0.0);
    L.reset();
    Tau.reset();
    Lt = std::make_shared<Tensor2d>("Lambda <iE|aF>", naoccB, navirA, navirB, navirA);
    Lt->sort(1243, Y, 1.0, 0.0);
    Y.reset();
    Lt->write(psio_, PSIF_DFOCC_AMPS);
    //Lt->print();
    Lt.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Tau (Q,IA) Intermediates //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Tau(Q,IA) = \sum{J,B} Tau(IJ,AB) b(Q,JB) + \sum{j,b} Tau(Ij,Ab) b(Q,jb)       (73)
    //\sum{J,B} Tau(IJ,AB) b(Q,JB)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS); 
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);
    T2.reset();
    X = std::make_shared<Tensor2d>("X (JB|IA)", naoccA, navirA, naoccA, navirA);
    X->sort(2413, Tau, 1.0, 0.0);
    Tau.reset();
    T = std::make_shared<Tensor2d>("Tau (Q|IA)", nQ, naoccA, navirA);
    T->gemm(false, false, bQiaA, X, 1.0, 0.0);
    X.reset();
    //+ \sum{j,b} Tau(Ij,Ab) b(Q,jb)
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    X = std::make_shared<Tensor2d>("X (jb|IA)", naoccB, navirB, naoccA, navirA);
    X->sort(2413, Tau, 1.0, 0.0);
    Tau.reset();
    T->gemm(false, false, bQiaB, X, 1.0, 1.0);
    X.reset();
    T->write(psio_, PSIF_DFOCC_AMPS);
    //T->print();
    T.reset();

    // Tau(Q,ia) = \sum{j,b} Tau(ij,ab) b(Q,jb) + \sum{J,B} Tau(Ji,Ba) b(Q,JB)       (74)
    //\sum{j,b} Tau(ij,ab) b(Q,jb)
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
    T2.reset();
    X = std::make_shared<Tensor2d>("X (jb|ia)", naoccB, navirB, naoccB, navirB);
    X->sort(2413, Tau, 1.0, 0.0);
    Tau.reset();
    T = std::make_shared<Tensor2d>("Tau (Q|ia)", nQ, naoccB, navirB);
    T->gemm(false, false, bQiaB, X, 1.0, 0.0);
    X.reset();
    //+ \sum{J,B} Tau(Ji,Ba) b(Q,JB)
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    X = std::make_shared<Tensor2d>("X (JB|ia)", naoccA, navirA, naoccB, navirB);
    X->sort(1324, Tau, 1.0, 0.0);
    Tau.reset();
    T->gemm(false, false, bQiaA, X, 1.0, 1.0);
    X.reset();
    T->write(psio_, PSIF_DFOCC_AMPS);
    //T->print();
    T.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Lambda_Q Intermediate /////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // L(Q) = \sum{M,E} L(M,E) b(Q,ME) + \sum{m,e} L(m,e) b(Q,me)       (76)
    L1c->gemv(false, bQiaA, l1A, 1.0, 0.0);
    L1c->gemv(false, bQiaB, l1B, 1.0, 1.0);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // L(Q,AI) Intermediates /////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // L(Q,AI) = \sum{M} L(M,A) b(Q,MI)        (78)
    T = std::make_shared<Tensor2d>("L1 (Q|AI)", nQ, navirA, naoccA);
    T->contract233(true, false, navirA, naoccA, l1A, bQijA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();
    // L(Q,ai) = \sum{m} L(m,a) b(Q,mi)        (79)
    T = std::make_shared<Tensor2d>("L1 (Q|ai)", nQ, navirB, naoccB);
    T->contract233(true, false, navirB, naoccB, l1B, bQijB, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gtilde(Q) Intermediate ////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gtilde(Q) = \sum{M,N} G(M,N) b(Q,MN) + \sum{m,n} G(m,n) b(Q,mn)       (81)
    gQt = std::make_shared<Tensor1d>("CCSD PDM G_Qt", nQ);
    gQt->gemv(false, bQijA, GijA, 1.0, 0.0);
    gQt->gemv(false, bQijB, GijB, 1.0, 1.0);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // G(Q,IJ) Intermediates /////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // G(Q,IJ) = \sum{M} G(I,M) b(Q,JM)        (83)
    T = std::make_shared<Tensor2d>("G (Q|IJ)", nQ, naoccA, naoccA);
    T->contract233(false, false, naoccA, naoccA, GijA, bQijA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();
    // G(Q,ij) = \sum{m} G(i,m) b(Q,jm)        (84)
    T = std::make_shared<Tensor2d>("G (Q|ij)", nQ, naoccB, naoccB);
    T->contract233(false, false, naoccB, naoccB, GijB, bQijB, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gtilde(Q,IA) Intermediates ////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gtilde(Q,IA) = \sum{M} Gtilde(M,I) b(Q,MA)        (86)
    T = std::make_shared<Tensor2d>("Gt (Q|IA)", nQ, naoccA, navirA);
    T->contract233(true, false, naoccA, navirA, GtijA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();
    // Gtilde(Q,ia) = \sum{m} Gtilde(m,i) b(Q,ma)        (87)
    T = std::make_shared<Tensor2d>("Gt (Q|ia)", nQ, naoccB, navirB);
    T->contract233(true, false, naoccB, navirB, GtijB, bQiaB, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();
 
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gtilde(Q,AI) Intermediates ////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gtilde(Q,AI) = \sum{E} Gtilde(A,E) b(Q,IE)        (89)
    T = std::make_shared<Tensor2d>("Gt (Q|AI)", nQ, navirA, naoccA);
    T->contract233(false, true, navirA, naoccA, GtabA, bQiaA, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();
    // Gtilde(Q,ai) = \sum{e} Gtilde(a,e) b(Q,ie)        (90)
    T = std::make_shared<Tensor2d>("Gt (Q|ai)", nQ, navirB, naoccB);
    T->contract233(false, true, navirB, naoccB, GtabB, bQiaB, 1.0, 0.0);
    T->write(psio_, PSIF_DFOCC_AMPS);
    T.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // V(Q,AB) and Vtilde(Q,AB) Intermediates ////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Alpha Block

    // V(Q,AB) = \sum{M,N} V(MA,NB) b(Q,MN) + \sum{m,n} V(mA,nB) b(Q,mn)       (92)
    //\sum{M,N} V(MA,NB) b(Q,MN)  1. terim
    V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (MN|AB)", naoccA, naoccA, navirA, navirA);
    X->sort(1324, V, 1.0, 0.0);
    V.reset();
    Vab = std::make_shared<Tensor2d>("V (Q|AB)", nQ, navirA, navirA);
    Vab->gemm(false, false, bQijA, X, 1.0, 0.0);
    //X.reset(); // Vtilde(Q,AB) nin 1.terimini hesaplarken kullaniyorum

    // Vtilde(Q,AB) = \sum{M,N} V(MA,NB) t(Q,NM) + \sum{m,n} V(mA,nB) t(Q,nm)   (95)
    //\sum{M,N} V(MA,NB) t(Q,NM)  1. terim
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T1 (Q|JI)", nQ, naoccA, naoccA);
    L->swap_3index_col(T);
    T.reset();
    Vtab = std::make_shared<Tensor2d>("Vt (Q|AB)", nQ, navirA, navirA);
    Vtab->gemm(false, false, L, X, 1.0, 0.0);
    X.reset();
    L.reset();

    // V(Q,AB) = \sum{M,N} V(MA,NB) b(Q,MN) + \sum{m,n} V(mA,nB) b(Q,mn)       (92)
    //+ \sum{m,n} V(mA,nB) b(Q,mn) 2. terim
    V = std::make_shared<Tensor2d>("V (iA|jB)", naoccB, navirA, naoccB, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (mn|AB)", naoccB, naoccB, navirA, navirA);
    X->sort(1324, V, 1.0, 0.0);
    V.reset();
    Vab->gemm(false, false, bQijB, X, 1.0, 1.0);
    //X.reset();  // Vtilde(Q,AB) nin 2.terimini hesaplarken kullaniyorum
    Vab->write(psio_, PSIF_DFOCC_AMPS);
    Vab.reset();

    // Vtilde(Q,AB) = \sum{M,N} V(MA,NB) t(Q,NM) + \sum{m,n} V(mA,nB) t(Q,nm)   (95)
    //+ \sum{m,n} V(mA,nB) t(Q,nm) 2. terim
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T1 (Q|ji)", nQ, naoccB, naoccB);
    L->swap_3index_col(T);
    T.reset();
    Vtab->gemm(false, false, L, X, 1.0, 1.0);
    Vtab->write(psio_, PSIF_DFOCC_AMPS);
    Vtab.reset();

// Beta Block

    // V(Q,ab) = \sum{m,n} V(ma,nb) b(Q,mn) + \sum{M,N} V(Ma,Nb) b(Q,MN)       (93)
    //\sum{m,n} V(ma,nb) b(Q,mn) 1.terim
    V = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (mn|ab)", naoccB, naoccB, navirB, navirB);
    X->sort(1324, V, 1.0, 0.0); 
    V.reset();
    Vab = std::make_shared<Tensor2d>("V (Q|ab)", nQ, navirB, navirB);
    Vab->gemm(false, false, bQijB, X, 1.0, 0.0);
    //X.reset(); // Vtilde (Q,ab) nin 1.terimini hesaplarken kullaniyorum

    // Vtilde(Q,ab) = \sum{m,n} V(ma,nb) t(Q,nm) + \sum{M,N} V(Ma,Nb) b(Q,nm)   (96)
    // \sum{m,n} V(ma,nb) t(Q,nm) // 1. terim
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T1 (Q|ji)", nQ, naoccB, naoccB);
    L->swap_3index_col(T);
    T.reset();
    Vtab = std::make_shared<Tensor2d>("Vt (Q|ab)", nQ, navirB, navirB);
    Vtab->gemm(false, false, L, X, 1.0, 0.0);
    X.reset();
    L.reset();

    // V(Q,ab) = \sum{m,n} V(ma,nb) b(Q,mn) + \sum{M,N} V(Ma,Nb) b(Q,MN)       (93)
    //+ \sum{M,N} V(Ma,Nb) b(Q,MN) 2. terim
    V = std::make_shared<Tensor2d>("V (Ia|Jb)", naoccA, navirB, naoccA, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (MN|ab)", naoccA, naoccA, navirB, navirB);
    X->sort(1324, V, 1.0, 0.0);
    V.reset();
    Vab->gemm(false, false, bQijA, X, 1.0, 1.0);
    //X.reset(); // Vtilde (Q,ab) nin 2.terimini hesaplarken kullaniyorum
    Vab->write(psio_, PSIF_DFOCC_AMPS);
    Vab.reset(); 

    // Vtilde(Q,ab) = \sum{m,n} V(ma,nb) t(Q,nm) + \sum{M,N} V(Ma,Nb) t(Q,nm)   (96)
    //+= \sum{M,N} V(Ma,Nb) t(Q,nm) // 2. terim
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("T1 (Q|JI)", nQ, naoccA, naoccA);
    L->swap_3index_col(T);
    T.reset();
    Vtab->gemm(false, false, L, X, 1.0, 1.0);
    L.reset();
    X.reset();
    Vtab->write(psio_, PSIF_DFOCC_AMPS);
    Vtab.reset();
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Nu(Q,IJ) Intermediates ////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Nu(Q,IJ) = \sum{M,E} V(IM,JE) b(Q,ME) + \sum{m,e} V(Im,Je) b(Q,me)       (98)
    V = std::make_shared<Tensor2d>("V <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ME|IJ)", naoccA, navirA, naoccA, naoccA);
    X->sort(2413, V, 1.0, 0.0);
    V.reset();
    Vij = std::make_shared<Tensor2d>("calV (Q|IJ)", nQ, naoccA, naoccA);
    Vij->gemm(false, false, bQiaA, X, 1.0, 0.0);
    X.reset();
    V = std::make_shared<Tensor2d>("V <Ij|Ka>", naoccA, naoccB, naoccA, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (me|IJ)", naoccB, navirB, naoccA, naoccA);
    X->sort(2413, V, 1.0, 0.0);
    V.reset();
    Vij->gemm(false, false, bQiaB, X, 1.0, 1.0);
    X.reset();
    Vij->write(psio_, PSIF_DFOCC_AMPS);
    Vij.reset();
    // Nu(Q,ij) = \sum{m,e} V(im,je) b(Q,me) + \sum{M,E} V(iM,jE) b(Q,ME)       (99)
    V = std::make_shared<Tensor2d>("V <ij|ka>", naoccB, naoccB, naoccB, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (me|ij)", naoccB, navirB, naoccB, naoccB);
    X->sort(2413, V, 1.0, 0.0);
    V.reset();
    Vij = std::make_shared<Tensor2d>("calV (Q|ij)", nQ, naoccB, naoccB);
    Vij->gemm(false, false, bQiaB, X, 1.0, 0.0);
    X.reset();
    V = std::make_shared<Tensor2d>("V <iJ|kA>", naoccB, naoccA, naoccB, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ME|ij)", naoccA, navirA, naoccB, naoccB);
    X->sort(2413, V, 1.0, 0.0);
    V.reset();
    Vij->gemm(false, false, bQiaA, X, 1.0, 1.0);
    X.reset();
    Vij->write(psio_, PSIF_DFOCC_AMPS);
    Vij.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Eta(Q,IJ) Intermediates ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Eta(Q,IJ) = \sum{M,E} Lambda(IM,JE) b(Q,ME) + \sum{m,e} Lambda(Im,Je) b(Q,me)       (101)
    L = std::make_shared<Tensor2d>("L <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ME|IJ)", naoccA, navirA, naoccA, naoccA);
    X->sort(2413, L, 1.0, 0.0);
    L.reset();
    Z = std::make_shared<Tensor2d>("Eta (Q|IJ)", nQ, naoccA, naoccA);
    Z->gemm(false, false, bQiaA, X, 1.0, 0.0);
    X.reset();
    L = std::make_shared<Tensor2d>("L <Ij|Ka>", naoccA, naoccB, naoccA, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (me|IJ)", naoccB, navirB, naoccA, naoccA);
    X->sort(2413, L, 1.0, 0.0);
    L.reset();
    Z->gemm(false, false, bQiaB, X, 1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();
    // Eta(Q,ij) = \sum{m,e} Lambda(im,je) b(Q,me) + \sum{M,E} Lambda(iM,jE) b(Q,ME)       (102)
    L = std::make_shared<Tensor2d>("L <ij|ka>", naoccB, naoccB, naoccB, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (me|ij)", naoccB, navirB, naoccB, naoccB);
    X->sort(2413, L, 1.0, 0.0);
    L.reset();
    Z = std::make_shared<Tensor2d>("Eta (Q|ij)", nQ, naoccB, naoccB);
    Z->gemm(false, false, bQiaB, X, 1.0, 0.0);    
    X.reset();
    L = std::make_shared<Tensor2d>("L <iJ|kA>", naoccB, naoccA, naoccB, navirA);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ME|ij)", naoccA, navirA, naoccB, naoccB);
    X->sort(2413, L, 1.0, 0.0);
    L.reset();
    Z->gemm(false, false, bQiaA, X, 1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Eta(Q,IA) Intermediates ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Eta(Q,IA) = \sum{M,N} Lambda(IM,NA) b(Q,MN) + \sum{m,n} Lambda(Im,nA) b(Q,mn)       (104)
    // Lambda(Ij,kA) = - Lambda(jI,kA)  ------->  Lambda(Im,nA) = - Lambda(mI,nA)
    // Eta(Q,IA) = \sum{M,N} Lambda(IM,NA) b(Q,MN) - \sum{m,n} Lambda(mI,nA)  b(Q,mn) 
    L = std::make_shared<Tensor2d>("L <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (MN|IA)", naoccA, naoccA, naoccA, navirA);
    X->sort(2314, L, 1.0, 0.0);
    L.reset();
    Z = std::make_shared<Tensor2d>("Eta (Q|IA)", nQ, naoccA, navirA);
    Z->gemm(false, false, bQijA, X, 1.0, 0.0);
    X.reset();
    L = std::make_shared<Tensor2d>("L <iJ|kA>", naoccB, naoccA, naoccB, navirA);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (mn|IA)", naoccB, naoccB, naoccA, navirA);
    X->sort(1324, L, 1.0, 0.0);
    L.reset();
    Z->gemm(false, false, bQijB, X, -1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();
    // Eta(Q,ia) = \sum{m,n} Lambda(im,na) b(Q,mn) + \sum{M,N} Lambda(iM,Na) b(Q,MN)       (105)
    // Lambda(iJ,Ka) = - Lambda(Ji,Ka)  -------> Lambda(iM,Na) = - Lambda(Mi,Na)
    // Eta(Q,ia) = \sum{m,n} Lambda(im,na) b(Q,mn) - \sum{M,N} Lambda(Mi,Na) b(Q,MN)
    L = std::make_shared<Tensor2d>("L <ij|ka>", naoccB, naoccB, naoccB, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (mn|ia)", naoccB, naoccB, naoccB, navirB);
    X->sort(2314, L, 1.0, 0.0);
    L.reset();
    Z = std::make_shared<Tensor2d>("Eta (Q|ia)", nQ, naoccB, navirB);
    Z->gemm(false, false, bQijB, X, 1.0, 0.0);    
    X.reset();
    L = std::make_shared<Tensor2d>("L <Ij|Ka>", naoccA, naoccB, naoccA, navirB);
    L->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (MN|ia)", naoccA, naoccA, naoccB, navirB);
    X->sort(1324, L, 1.0, 0.0);
    L.reset();
    Z->gemm(false, false, bQijA, X, -1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Etatilde(Q,IA) Intermediates //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Etatilde(Q,IA) = \sum{E,F} Lambdatilde(IE,AF) b(Q,EF) + \sum{e,f} Lambdatilde(Ie,Af) b(Q,ef)       (107)
    Lt = std::make_shared<Tensor2d>("Lambda <IE|AF>", naoccA, navirA, navirA, navirA);
    Lt->read(psio_, PSIF_DFOCC_AMPS);    
    X = std::make_shared<Tensor2d>("X (EF|IA)", navirA, navirA, naoccA, navirA);
    X->sort(2413, Lt, 1.0, 0.0);
    Lt.reset();
    Z = std::make_shared<Tensor2d>("Eta2 (Q|IA)", nQ, naoccA, navirA);
    Z->gemm(false, false, bQabA, X, 1.0, 0.0);
    X.reset();
    Lt = std::make_shared<Tensor2d>("Lambda <Ie|Af>", naoccA, navirB, navirA, navirB);
    Lt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ef|IA)", navirB, navirB, naoccA, navirA);
    X->sort(2413, Lt, 1.0, 0.0);
    Lt.reset();
    Z->gemm(false, false, bQabB, X, 1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();
    // Etatilde(Q,ia) = \sum{e,f} Lambdatilde(ie,af) b(Q,ef) + \sum{E,F} Lambdatilde(iE,aF) b(Q,EF)       (108)
    Lt = std::make_shared<Tensor2d>("Lambda <ie|af>", naoccB, navirB, navirB, navirB);
    Lt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ef|ia)", navirB, navirB, naoccB, navirB);
    X->sort(2413, Lt, 1.0, 0.0);
    Lt.reset();
    Z = std::make_shared<Tensor2d>("Eta2 (Q|ia)", nQ, naoccB, navirB);
    Z->gemm(false, false, bQabB, X, 1.0, 0.0);
    X.reset();
    Lt = std::make_shared<Tensor2d>("Lambda <iE|aF>", naoccB, navirA, navirB, navirA);
    Lt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (EF|ia)", navirA, navirA, naoccB, navirB);
    X->sort(2413, Lt, 1.0, 0.0);
    Lt.reset();
    Z->gemm(false, false, bQabA, X, 1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Etatilde(Q,AB) Intermediates //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Etatilde(Q,AB) = \sum{M,E} Lambdatilde(MA,EB) b(Q,ME) + \sum{m,e} Lambdatilde(mA,eB) b(Q,me)       (110)
    Lt = std::make_shared<Tensor2d>("Lambda <IE|AF>", naoccA, navirA, navirA, navirA);
    Lt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ME|AB)", naoccA, navirA, navirA, navirA);
    X->sort(1324, Lt, 1.0, 0.0);
    Lt.reset();
    Z = std::make_shared<Tensor2d>("Eta2 (Q|AB)", nQ, navirA, navirA);
    Z->gemm(false, false, bQiaA, X, 1.0, 0.0);
    X.reset();
    Lt = std::make_shared<Tensor2d>("Lambda <iE|aF>", naoccB, navirA, navirB, navirA);
    Lt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (me|AB)", naoccB, navirB, navirA, navirA);
    X->sort(1324, Lt, 1.0, 0.0);
    Lt.reset();
    Z->gemm(false, false, bQiaB, X, 1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();
    // Etatilde(Q,ab) = \sum{m,e} Lambdatilde(ma,eb) b(Q,me) + \sum{M,E} Lambdatilde(Ma,Eb) b(Q,ME)       (111)
    Lt = std::make_shared<Tensor2d>("Lambda <ie|af>", naoccB, navirB, navirB, navirB);
    Lt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (me|ab)", naoccB, navirB, navirB, navirB);
    X->sort(1324, Lt, 1.0, 0.0);
    Lt.reset();
    Z = std::make_shared<Tensor2d>("Eta2 (Q|ab)", nQ, navirB, navirB);
    Z->gemm(false, false, bQiaB, X, 1.0, 0.0);
    X.reset();
    Lt = std::make_shared<Tensor2d>("Lambda <Ie|Af>", naoccA, navirB, navirA, navirB);
    Lt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ME|ab)", naoccA, navirA, navirB, navirB);
    X->sort(1324, Lt, 1.0, 0.0);
    Lt.reset();
    Z->gemm(false, false, bQiaA, X, 1.0, 1.0);
    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
}  // end ccsd_pdm_3index_intr()


void DFOCC::uccsd_pdm_yQia() { 

    SharedTensor2d T, T2, L2, Tau, X, Y, Z, V, Vt, U, L, Lt, Yt;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Y (IA,JB)  Intermediates //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // AAAA Block
    // Y(IA,JB) = \sum{M,E} [t(MI,AE) + 2 {t(M,A) * t(I,E)}] * Vtilde(JE,MB) -  \sum{m,e} t(Im,Ae) * Vtilde(Je,mB)       (55)
    X = std::make_shared<Tensor2d>("X (MA|IE)", naoccA, navirA, naoccA, navirA);
    X->dirprd224(t1A, t1A, 2.0, 0.0);
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X->sort(1324, T2, 1.0, 1.0);
    T2.reset();
    U = std::make_shared<Tensor2d>("X (IA|ME)", naoccA, navirA, naoccA, navirA);
    U->sort(3214, X, 1.0, 0.0);
    X.reset();
    Vt = std::make_shared<Tensor2d>("Vt (IA|JB)", naoccA, navirA, naoccA, navirA);
    Vt->read(psio_, PSIF_DFOCC_AMPS);    
    V = std::make_shared<Tensor2d>("Vtx (ME|JB)", naoccA, navirA, naoccA, navirA);
    V->sort(3214, Vt, 1.0, 0.0);
    Vt.reset();
    Y = std::make_shared<Tensor2d>("Y (IA|JB)", naoccA, navirA, naoccA, navirA);
    Y->gemm(false, false, U, V, 1.0, 0.0);
    U.reset();
    V.reset();

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (IA|me)", naoccA, navirA, naoccB, navirB);
    X->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    Vt = std::make_shared<Tensor2d>("Vt (Ia|jB)", naoccA, navirB, naoccB, navirA);
    Vt->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("Vtx (me|JB)", naoccB, navirB, naoccA, navirA);
    V->sort(3214, Vt, 1.0, 0.0);
    Vt.reset();
    Y->gemm(false, false, X, V, -1.0, 1.0);
    V.reset();
    X.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
    Y.reset();
 
    // BBBB Block
    // Y(ia,jb) = \sum{m,e} [t(mi,ae) + 2 {t(m,a) * t(i,e)}] * Vtilde(je,mb) -  \sum{M,E} t(Mi,Ea) * Vtilde(jE,Mb)       (56)
    X = std::make_shared<Tensor2d>("X (ma|ie)", naoccB, navirB, naoccB, navirB);
    X->dirprd224(t1B, t1B, 2.0, 0.0);
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X->sort(1324, T2, 1.0, 1.0);
    T2.reset();
    U = std::make_shared<Tensor2d>("X (ia|me)", naoccB, navirB, naoccB, navirB);
    U->sort(3214, X, 1.0, 0.0);
    X.reset();
    Vt = std::make_shared<Tensor2d>("Vt (ia|jb)", naoccB, navirB, naoccB, navirB);
    Vt->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("Vtx (me|jb)", naoccB, navirB, naoccB, navirB);
    V->sort(3214, Vt, 1.0, 0.0);
    Vt.reset();
    Y = std::make_shared<Tensor2d>("Y (ia|jb)", naoccB, navirB, naoccB, navirB);
    Y->gemm(false, false, U, V, 1.0, 0.0);
    U.reset();
    V.reset();

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ia|ME)", naoccB, navirB, naoccA, navirA);
    X->sort(2413, T2, 1.0, 0.0);
    T2.reset();
    Vt = std::make_shared<Tensor2d>("Vt (iA|Jb)", naoccB, navirA, naoccA, navirB);
    Vt->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("Vtx (ME|jb)", naoccA, navirA, naoccB, navirB);
    V->sort(3214, Vt, 1.0, 0.0);
    Vt.reset();
    Y->gemm(false, false, X, V, -1.0, 1.0);
    V.reset();
    X.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
    Y.reset();

    // AABB Block
    // Y_{IAjb} = \sum_{M,E}(t_{MI,AE} + 2 t_{M,A} t_{I,E}) \widetilde{V}_{jEMb} - \sum_{m,e} t_{Im}^{Ae} \widetilde{V}_{jemb}
    X = std::make_shared<Tensor2d>("X (MA|IE)", naoccA, navirA, naoccA, navirA);
    X->dirprd224(t1A, t1A, 2.0, 0.0);
//t1A->print();
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X->sort(1324, T2, 1.0, 1.0);
    T2.reset();
    U = std::make_shared<Tensor2d>("X (IA|ME)", naoccA, navirA, naoccA, navirA);
    U->sort(3214, X, 1.0, 0.0);
    X.reset();
    Vt = std::make_shared<Tensor2d>("Vt (iA|Jb)", naoccB, navirA, naoccA, navirB);
    Vt->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("Vtx (ME|jb)", naoccA, navirA, naoccB, navirB);
    V->sort(3214, Vt, 1.0, 0.0);
//V->print();
    Vt.reset();
    Y = std::make_shared<Tensor2d>("Y (IA|jb)", naoccA, navirA, naoccB, navirB);
    Y->gemm(false, false, U, V, 1.0, 0.0);
    U.reset();
    V.reset();

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
/*
SharedTensor2d A = std::make_shared<Tensor2d>("A (IA|jb)", naoccA, navirA, naoccB, navirB);
A->sort(1324, T2, 1.0, 0.0);
A->print();
A.reset();
*/
    X = std::make_shared<Tensor2d>("X (IA|me)", naoccA, navirA, naoccB, navirB);
    X->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    Vt = std::make_shared<Tensor2d>("Vt (ia|jb)", naoccB, navirB, naoccB, navirB);
    Vt->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("Vtx (me|jb)", naoccB, navirB, naoccB, navirB);
    V->sort(3214, Vt, 1.0, 0.0);
//V->print();
    Vt.reset();
    Y->gemm(false, false, X, V, -1.0, 1.0);
    V.reset();
    X.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
//Y->print();
    Y.reset();

    // BBAA Block
    // Y_{iaJB} = \sum_{m,e} (t_{mi,ae} + 2 t_{m,a} t_{i,e}) \widetilde{V}_{JemB} - \sum_{M,E}(t_{Mi}^{Ea} \widetilde{V}_{JEMB}
    X = std::make_shared<Tensor2d>("X (ma|ie)", naoccB, navirB, naoccB, navirB);
    X->dirprd224(t1B, t1B, 2.0, 0.0);
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X->sort(1324, T2, 1.0, 1.0);
    T2.reset();
    U = std::make_shared<Tensor2d>("X (ia|me)", naoccB, navirB, naoccB, navirB);
    U->sort(3214, X, 1.0, 0.0);
    X.reset();
    Vt = std::make_shared<Tensor2d>("Vt (Ia|jB)", naoccA, navirB, naoccB, navirA);
    Vt->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("Vtx (me|JB)", naoccB, navirB, naoccA, navirA);
    V->sort(3214, Vt, 1.0, 0.0);
    Vt.reset();
    Y = std::make_shared<Tensor2d>("Y (ia|JB)", naoccB, navirB, naoccA, navirA);
    Y->gemm(false, false, U, V, 1.0, 0.0);
    U.reset();
    V.reset();

    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ia|ME)", naoccB, navirB, naoccA, navirA);
    X->sort(2413, T2, 1.0, 0.0);
    T2.reset();
    Vt = std::make_shared<Tensor2d>("Vt (IA|JB)", naoccA, navirA, naoccA, navirA);
    Vt->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("Vtx (ME|JB)", naoccA, navirA, naoccA, navirA);
    V->sort(3214, Vt, 1.0, 0.0);
    Vt.reset();
    Y->gemm(false, false, X, V, -1.0, 1.0);
    V.reset();
    X.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
//Y->print();
    Y.reset();

    // ABBA Block
    // Y_{IajB} = \sum_{m,E} (t_{Im,Ea} + 2 t_{m,a} t_{I,E}) \widetilde{V}_{jEmB} 
    X = std::make_shared<Tensor2d>("X (ma|IE)", naoccB, navirB, naoccA, navirA);
    X->dirprd224(t1B, t1A, 2.0, 0.0);
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    X->sort(2413, T2, 1.0, 1.0);
    T2.reset();
    U = std::make_shared<Tensor2d>("X (Ia|mE)", naoccA, navirB, naoccB, navirA);
    U->sort(3214, X, 1.0, 0.0);
    X.reset(); 
    Vt = std::make_shared<Tensor2d>("Vt (iA|jB)", naoccB, navirA, naoccB, navirA);
    Vt->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("Vtx (mE|jB)", naoccB, navirA, naoccB, navirA);
    V->sort(3214, Vt, 1.0, 0.0);
    Vt.reset();
    Y = std::make_shared<Tensor2d>("Y (Ia|jB)", naoccA, navirB, naoccB, navirA);
    Y->gemm(false, false, U, V, 1.0, 0.0);
    U.reset();
    V.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
    Y.reset();

    // BAAB Block
    // Y_{iAJb} = \sum_{M,e} (t_{Mi,Ae} + 2 t_{M,A} t_{i,e}) \widetilde{V}_{JeMb}
    X = std::make_shared<Tensor2d>("X (MA|ie)", naoccA, navirA, naoccB, navirB);
    X->dirprd224(t1A, t1B, 2.0, 0.0);
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    X->sort(1324, T2, 1.0, 1.0);
    T2.reset();
    U = std::make_shared<Tensor2d>("X (iA|Me)", naoccB, navirA, naoccA, navirB);
    U->sort(3214, X, 1.0, 0.0);
    X.reset();
    Vt = std::make_shared<Tensor2d>("Vt (Ia|Jb)", naoccA, navirB, naoccA, navirB);
    Vt->read(psio_, PSIF_DFOCC_AMPS);
    V = std::make_shared<Tensor2d>("Vtx (Me|Jb)", naoccA, navirB, naoccA, navirB);
    V->sort(3214, Vt, 1.0, 0.0);
    Vt.reset();
    Y = std::make_shared<Tensor2d>("Y (iA|Jb)", naoccB, navirA, naoccA, navirB);
    Y->gemm(false, false, U, V, 1.0, 0.0);
    U.reset();
    V.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
    Y.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Y (IJ,AB)  Intermediates //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Y(IJ,AB) = 0.5 * \sum{M,N} Tau(MN,AB) * V(IJ,MN)       (60)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);
    T2.reset();
    V = std::make_shared<Tensor2d>("V <IJ|KL>", naoccA, naoccA, naoccA, naoccA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("Y <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Y->gemm(false, false, V, Tau, 0.5, 0.0);
    Tau.reset();
    V.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
    Y.reset();
 
    // Y(ij,ab) = 0.5 * \sum{m,n} Tau(mn,ab) * V(ij,mn)       (61)
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
    T2.reset();
    V = std::make_shared<Tensor2d>("V <ij|kl>", naoccB, naoccB, naoccB, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("Y <ij|ab>", naoccB, naoccB, navirB, navirB);
    Y->gemm(false, false, V, Tau, 0.5, 0.0);
    Tau.reset();
    V.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
    Y.reset();

    // Y(Ij,Ab) = \sum{M,n} Tau(Mn,Ab) * V(Ij,Mn)        (62)
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    V = std::make_shared<Tensor2d>("V <Ij|Kl>", naoccA, naoccB, naoccA, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("Y <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Y->gemm(false, false, V, Tau, 1.0, 0.0);
    Tau.reset();
    V.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
/*
SharedTensor2d A = std::make_shared<Tensor2d>("Y <IJ|AB>", naoccA, naoccA, navirA, navirA);
A->read(psio_, PSIF_DFOCC_AMPS);
A->axpy(Y, 1.0);
A->print();
A.reset();
*/
    Y.reset();

    // Y(iJ,aB) = \sum{m,N} Tau(Nm,Ba) * V(iJ,mN) - 0.5 * \sum{M,n} Tau(Mn,Ba) * V(iJ,Mn)       (63)   // hatali sonucta kontrol et.
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    X = std::make_shared<Tensor2d>("X <mN|aB>", naoccB, naoccA, navirB, navirA);
    X->sort(2143, Tau, 1.0, 0.0);
    Tau.reset();
    V = std::make_shared<Tensor2d>("V <Ij|Kl>", naoccA, naoccB, naoccA, naoccB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    U = std::make_shared<Tensor2d>("X <iJ|mN>", naoccB, naoccA, naoccB, naoccA);
    U->sort(2143, V, 1.0, 0.0);
    V.reset();
    Y = std::make_shared<Tensor2d>("Y <iJ|aB>", naoccB, naoccA, navirB, navirA);
    Y->gemm(false, false, U, X, 1.0, 0.0);
    X.reset();
    U.reset();
    Y->write(psio_, PSIF_DFOCC_AMPS);
    Y.reset();

/*
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // s (I,A)  Intermediates //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // s(I,A) = \sum{M} t(M,A) \sum{E} t(I,E) * Lambda(M,E)       (70)
    SharedTensor2d Sij = std::make_shared<Tensor2d>("S <I|J>", naoccA, naoccA);
    SharedTensor2d s1A = std::make_shared<Tensor2d>("S <I|A>", naoccA, navirA);
    Sij->gemm(false, true, t1A, l1A, 1.0, 0.0);
    s1A->gemm(false, false, Sij, t1A, 1.0, 0.0);      // Ytilde <IJ|AB> teriminin sonunda resetliyorum.
    Sij.reset();
    // s(i,a) = \sum{m} t(m,a) \sum{e} t(i,e) * Lambda(m,e)       (71)
    Sij = std::make_shared<Tensor2d>("S <i|j>", naoccB, naoccB);
    SharedTensor2d s1B = std::make_shared<Tensor2d>("S <i|a>", naoccB, navirB);
    Sij->gemm(false, true, t1B, l1B, 1.0, 0.0);
    s1B->gemm(false, false, Sij, t1B, 1.0, 0.0);      // Ytilde <ij|ab> teriminin sonunda resetliyorum.
    Sij.reset();
*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Ytilde (IJ,AB)  Intermediates /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Ytilde(IJ,AB) = 0.5 * [Y(IJ,AB) - Y(IA,JB) + Y(JA,IB) + Y(IB,JA) - Y(JB,IA)] 
                    // 1.5 * [t(I,A)* s(J,B) - t(J,A)* s(I,B) - t(I,B)* s(J,A) + t(J,B)* s(I,A)] + V(IB,JA) + V(JA,IB)    (65) 

    Y = std::make_shared<Tensor2d>("Y <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt = std::make_shared<Tensor2d>("Y2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Yt->axpy(Y, 0.5);
    Y.reset();
    U = std::make_shared<Tensor2d>("Y (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(1324, U, -0.5, 1.0);
    Yt->sort(3124, U, 0.5, 1.0);
    Yt->sort(1342, U, 0.5, 1.0);
    Yt->sort(3142, U, -0.5, 1.0);
    U.reset();

    // s(I,A) = \sum{M} t(M,A) \sum{E} t(I,E) * Lambda(M,E)       (70)
    SharedTensor2d Sij = std::make_shared<Tensor2d>("S <I|J>", naoccA, naoccA);
    SharedTensor2d s1A = std::make_shared<Tensor2d>("S <I|A>", naoccA, navirA);
    Sij->gemm(false, true, t1A, l1A, 1.0, 0.0);
    s1A->gemm(false, false, Sij, t1A, 1.0, 0.0);      // Ytilde <IJ|AB> teriminin sonunda resetliyorum.
    Sij.reset();

    // Ytilde(IJ,AB) += 1.5 * [t(I,A)* s(J,B) - t(J,A)* s(I,B) - t(I,B)* s(J,A) + t(J,B)* s(I,A)] 
    for (int i = 0; i < naoccA; ++i) {
        for (int j = 0; j < naoccA; ++j) {
             int ij = (i * naoccA) + j;
            for (int a = 0; a < navirA; ++a) {
                for (int b = 0; b < navirA; ++b) {
                    int ab = (a * navirA) + b;
                    double value = (t1A->get(i, a) * s1A->get(j, b)) - (t1A->get(j, a) * s1A->get(i, b)) 
							- (t1A->get(i, b) * s1A->get(j, a)) + (t1A->get(j, b) * s1A->get(i, a));
                    Yt->add(ij, ab, 1.5 * value);
                }
            }
        }
    }
    
    V = std::make_shared<Tensor2d>("V (IA|JB)", naoccA, navirA, naoccA, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(1342, V, 1.0, 1.0);
    Yt->sort(3124, V, 1.0, 1.0);
    V.reset();
    Yt->write(psio_, PSIF_DFOCC_AMPS);
//Yt->print();
    Yt.reset();

    // Ytilde(ij,ab) = 0.5 * [Y(ij,ab) - Y(ia,jb) + Y(ja,ib) + Y(ib,ja) - Y(jb,ia)] 
                    // 1.5 * [t(i,a)* s(j,b) - t(j,a)* s(i,b) - t(i,b)* s(j,a) + t(j,b)* s(i,a)] + V(ib,ja) + V(ja,ib)    (66) 
    Y = std::make_shared<Tensor2d>("Y <ij|ab>", naoccB, naoccB, navirB, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt = std::make_shared<Tensor2d>("Y2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    Yt->axpy(Y, 0.5);
    Y.reset();
    U = std::make_shared<Tensor2d>("Y (ia|jb)", naoccB, navirB, naoccB, navirB);
    U->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(1324, U, -0.5, 1.0);
    Yt->sort(3124, U, 0.5, 1.0);
    Yt->sort(1342, U, 0.5, 1.0);
    Yt->sort(3142, U, -0.5, 1.0);
    U.reset();

    // s(i,a) = \sum{m} t(m,a) \sum{e} t(i,e) * Lambda(m,e)       (71)
    Sij = std::make_shared<Tensor2d>("S <i|j>", naoccB, naoccB);
    SharedTensor2d s1B = std::make_shared<Tensor2d>("S <i|a>", naoccB, navirB);
    Sij->gemm(false, true, t1B, l1B, 1.0, 0.0);
    s1B->gemm(false, false, Sij, t1B, 1.0, 0.0);      // Ytilde <ij|ab> teriminin sonunda resetliyorum.
    Sij.reset();
 
    for (int i = 0; i < naoccB; ++i) {
        for (int j = 0; j < naoccB; ++j) {
             int ij = (i * naoccB) + j; 
            for (int a = 0; a < navirB; ++a) {
                for (int b = 0; b < navirB; ++b) {
                    int ab = (a * navirB) + b;
                    double value = (t1B->get(i, a) * s1B->get(j, b)) - (t1B->get(j, a) * s1B->get(i, b)) 
							- (t1B->get(i, b) * s1B->get(j, a)) + (t1B->get(j, b) * s1B->get(i, a));
                    Yt->add(ij, ab, 1.5 * value);
                }
            }
        }
    }

    V  = std::make_shared<Tensor2d>("V (ia|jb)", naoccB, navirB, naoccB, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(1342, V, 1.0, 1.0);
    Yt->sort(3124, V, 1.0, 1.0);
    V.reset();
    Yt->write(psio_, PSIF_DFOCC_AMPS);
    Yt.reset();

    // Ytilde(Ij,Ab) = 0.5 * [Y(Ij,Ab) - Y(IA,jb) + Y(jA,Ib) + Y(Ib,jA) - Y(jb,IA)] 
                    // 1.5 * [t(I,A)* s(j,b) + t(j,b)* s(I,A)] + V(Ib,jA) + V(jA,Ib)     (67) 
    Y = std::make_shared<Tensor2d>("Y <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
//Y->print();
    Yt = std::make_shared<Tensor2d>("Y2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Yt->axpy(Y, 0.5);
    Y.reset();
    Y = std::make_shared<Tensor2d>("Y (IA|jb)", naoccA, navirA, naoccB, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
//Y->print();
    Yt->sort(1324, Y, -0.5, 1.0);
    Y.reset();
    Y = std::make_shared<Tensor2d>("Y (iA|Jb)", naoccB, navirA, naoccA, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(3124, Y, 0.5, 1.0);
    Y.reset();
//Yt->print();
    Y = std::make_shared<Tensor2d>("Y (Ia|jB)", naoccA, navirB, naoccB, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(1342, Y, 0.5, 1.0);
    Y.reset();
    Y = std::make_shared<Tensor2d>("Y (ia|JB)", naoccB, navirB, naoccA, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(3142, Y, -0.5, 1.0);
    Y.reset();
//Yt->print();
/* 
    Y = std::make_shared<Tensor2d>("Y <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    SharedTensor2d Yt2 = std::make_shared<Tensor2d>("Y2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Yt2->axpy(Y, 0.5);
    Y.reset();
    U = std::make_shared<Tensor2d>("Y (IA|JB)", naoccA, navirA, naoccA, navirA);
    U->read(psio_, PSIF_DFOCC_AMPS);
    Yt2->sort(1324, U, -0.5, 1.0);
    Yt2->sort(3124, U, 0.5, 1.0);
    Yt2->sort(1342, U, 0.5, 1.0);
    Yt2->sort(3142, U, -0.5, 1.0);
    U.reset();
Yt2->axpy(Yt, 1.0);
Yt2->print();
*/

    // 1.5 * [t(I,A)* s(j,b) + t(j,b)* s(I,A)] 
    for (int i = 0; i < naoccA; ++i) {
        for (int j = 0; j < naoccB; ++j) {
             int ij = (i * naoccB) + j;
            for (int a = 0; a < navirA; ++a) {
                for (int b = 0; b < navirB; ++b) {
                    int ab = (a * navirB) + b;
                    double value = (t1A->get(i, a) * s1B->get(j, b)) + (t1B->get(j, b) * s1A->get(i, a));
                    Yt->add(ij, ab, 1.5 * value);
                }
            }
        }
    }

    V = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(1342, V, 1.0, 1.0);
    V.reset();
    V = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(3124, V, 1.0, 1.0);
    V.reset();
    Yt->write(psio_, PSIF_DFOCC_AMPS);
//Yt->print();
    Yt.reset();    

    // Ytilde(iJ,aB) = 0.5 * [Y(iJ,aB) - Y(ia,JB) + Y(Ja,iB) + Y(iB,Ja) - Y(JB,ia)] 
                    // 1.5 * [t(i,a)* s(J,B) - t(J,B)* s(i,a)] + V(iB,Ja) + V(Ja,iB)    (65) 
    Y = std::make_shared<Tensor2d>("Y <iJ|aB>", naoccB, naoccA, navirB, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt = std::make_shared<Tensor2d>("Y2 <iJ|aB>", naoccB, naoccA, navirB, navirA);
    Yt->axpy(Y, 0.5);
    Y.reset();

    Y = std::make_shared<Tensor2d>("Y (ia|JB)", naoccB, navirB, naoccA, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(1324, Y, -0.5, 1.0);
    Y.reset();

    Y = std::make_shared<Tensor2d>("Y (Ia|jB)", naoccA, navirB, naoccB, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(3124, Y, 0.5, 1.0);
    Y.reset();

    Y = std::make_shared<Tensor2d>("Y (iA|Jb)", naoccB, navirA, naoccA, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(1342, Y, 0.5, 1.0);
    Y.reset();

    Y = std::make_shared<Tensor2d>("Y (IA|jb)", naoccA, navirA, naoccB, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(3142, Y, -0.5, 1.0);
    Y.reset();


    // 1.5 * [t(i,a)* s(J,B) - t(J,B)* s(i,a)] 
    for (int i = 0; i < naoccB; ++i) {
        for (int j = 0; j < naoccA; ++j) {
             int ij = (i * naoccA) + j;
            for (int a = 0; a < navirB; ++a) {
                for (int b = 0; b < navirA; ++b) {
                    int ab = (a * navirA) + b;
                    double value = (t1B->get(i, a) * s1A->get(j, b)) + (t1A->get(j, b) * s1B->get(i, a));
                    Yt->add(ij, ab, 1.5 * value);
                }
            }
        }
    }
    s1A.reset();
    s1B.reset();

    V = std::make_shared<Tensor2d>("V (iA|Jb)", naoccB, navirA, naoccA, navirB);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(1342, V, 1.0, 1.0);
    V.reset();
    V = std::make_shared<Tensor2d>("V (Ia|jB)", naoccA, navirB, naoccB, navirA);
    V->read(psio_, PSIF_DFOCC_AMPS);
    Yt->sort(3124, V, 1.0, 1.0);
    V.reset();
    Yt->write(psio_, PSIF_DFOCC_AMPS);
//Yt->print();
    Yt.reset();

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // y(Q,IA) Intermediates /////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // y(Q,IA) = \sum{M,E} Ytilde(IM,AE) b(Q,ME) + \sum{m,e} Ytilde(Im,Ae) b(Q,me)       (113)
    Yt = std::make_shared<Tensor2d>("Y2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Yt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ME|IA)", naoccA, navirA, naoccA, navirA);
    X->sort(2413, Yt, 1.0, 0.0);
    Yt.reset();
    Z = std::make_shared<Tensor2d>("Y (Q|IA)", nQ, naoccA, navirA);
    Z->gemm(false, false, bQiaA, X, 1.0, 0.0);
    X.reset();
    Yt = std::make_shared<Tensor2d>("Y2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    Yt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (me|IA)", naoccB, navirB, naoccA, navirA);
    X->sort(2413, Yt, 1.0, 0.0);
    Yt.reset();
    Z->gemm(false, false, bQiaB, X, 1.0, 1.0);

    /* For RHF correction of Y2 <IJ|AB>
    Y = std::make_shared<Tensor2d>("Y2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    SharedTensor2d A = std::make_shared<Tensor2d>("X (ME|IA)", naoccA, navirA, naoccA, navirA);
    A->sort(2413, Y, 1.0, 0.0);
    Y.reset();
    A->axpy(X, 1.0);
    A->print();
    A.reset();
    */

    X.reset();
    Z->write(psio_, PSIF_DFOCC_AMPS);
//Z->print();
    Z.reset();

    // y(Q,ia) = \sum{m,e} Ytilde(im,ae) b(Q,me) + \sum{M,E} Ytilde(iM,aE) b(Q,ME)       (114)
    Yt = std::make_shared<Tensor2d>("Y2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    Yt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (me|ia)", naoccB, navirB, naoccB, navirB);
    X->sort(2413, Yt, 1.0, 0.0);
    Yt.reset();
    Z = std::make_shared<Tensor2d>("Y (Q|ia)", nQ, naoccB, navirB);
    Z->gemm(false, false, bQiaB, X, 1.0, 0.0);
    X.reset();
    Yt = std::make_shared<Tensor2d>("Y2 <iJ|aB>", naoccB, naoccA, navirB, navirA);
    Yt->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (ME|ia)", naoccA, navirA, naoccB, navirB);
    X->sort(2413, Yt, 1.0, 0.0);
    Yt.reset();
    Z->gemm(false, false, bQiaA, X, 1.0, 1.0);

    /* For RHF correction of Y2 <IJ|AB>
    Y = std::make_shared<Tensor2d>("Y2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    SharedTensor2d A = std::make_shared<Tensor2d>("X (me|ia)", naoccB, navirB, naoccB, navirB);
    A->sort(2413, Y, 1.0, 0.0);
    Y.reset();
    A->axpy(X, 1.0);
    A->print();
    A.reset();
    */

    X.reset();
//Z->print();
    Z->write(psio_, PSIF_DFOCC_AMPS);
} // end ccsd_pdm_yQia

}  // namespace dfoccwave
}  // namespace psi



