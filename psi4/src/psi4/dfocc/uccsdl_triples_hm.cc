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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <ctime>
#include "dfocc.h"
#include "defines.h"
#include "psi4/libqt/qt.h"

namespace psi {
namespace dfoccwave {

/*
 * You can jump anywhere in file. Just search N<number>
 * 0  : this list (N0)
 * 1  : beginning of AAA
 * 2  : main loop of AAA
 * 3  : energy of AAA
 * 4  : beginning of BBB
 * 5  : main loop of BBB
 * 6  : energy of BBB
 * 7  : beginning of AAB
 * 8  : main loop of AAB
 * 9  : energy of AAB
 * 10 : beginning of ABB
 * 11 : main loop of ABB
 * 12 : energy of ABB
 */

void DFOCC::uccsdl_triples_hm()
{
    pair_index();

    outfile->Printf("\tUsing high-memory disk algorithm...\n\n");

    timer_on("CCSD(AT)-HM-AAA");
    //==================================================
    //======================= AAA ======================
    //==================================================
    // N1 : beginning of AAA

    // form <IA||BC>
    // 'c' letter mean compact
    SharedTensor2d Jc_I_ABC = std::make_shared<Tensor2d>("J[I] (A|B>=C)", navirA, ntri_abAA);
    //SharedTensor2d J_I_ABC = std::make_shared<Tensor2d>("J[I] (A|BC)", navirA, navirA, navirA);
    SharedTensor2d G_I_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);
    SharedTensor2d bQabA_c = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    bQabA_c->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d biaQA = std::make_shared<Tensor2d>("B (IA|Q)", naoccA * navirA, nQ);
    biaQA->trans(bQiaA);

    for (long int i = 0; i < naoccA; i++) {
        bool flag = (i == 0) ? false : true;
        // [I](A|B>=C) = [I](A|Q) * (Q|A>=B)
        Jc_I_ABC->contract(false, false, navirA, ntri_abAA, nQ, biaQA, bQabA_c, i * navirA * nQ, 0, 1.0, 0.0);

        // expand [I](A|B>=C) to [I](A|BC)
        // [I]<B||AC> = [I]<B|AC> - [I]<B|CA>
        //            = [I](A|BC) - [I](C|BA)
        for (long int a = 0; a < navirA; a++) {
            for (long int b = 0; b < navirA; b++) {
                long int ab = ab_idxAA->get(a, b);
                long int ba = ab_idxAA->get(b, a);
                for (int c = 0; c < navirA; c++) {
                    long int ac = ab_idxAA->get(a, c);
                    long int ca = ab_idxAA->get(c, a);
                    long int bc = ab_idxAA->get(b, c);
                    //double val = Jc_I_ABC->get(a, index2(b, c)); // for only expand
                    //J_I_ABC->set(a, bc, val); // for only expand
                    //J_I_ABC->set(b, ac, val); // for only sort
                    double val = Jc_I_ABC->get(b, index2(a, c)) - Jc_I_ABC->get(c, index2(b, a)); // for expand, sort and subtract (anti-symmetrization)
                    G_I_ABC->set(a, bc, val); // for expand, sort and subtract (anti-symmetrization)
                } // c
            } // b
        } // a

        // write [I]<A||BC> to disk
        G_I_ABC->mywrite(psio_, PSIF_DFOCC_IABC_AAAA, flag);
    } // i

    Jc_I_ABC.reset();
    G_I_ABC.reset();
    bQabA_c.reset();
    //J_I_ABC.reset();

    // form <IJ||AK>
    // (IA|JK) = bQiaA * bQijA
    // <IJ||AK> = <AK||IJ> = <AK|IJ> - <AK|JI>
    //                     = (AI|KJ) - (AJ|KI)
    //                     = (IA|JK) - (JA|IK)
    // (JA|IK) = sort 3214 (IA|JK)
    SharedTensor2d J_IAJK = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JK)", naoccA, navirA, naoccA, naoccA);
    J_IAJK->gemm(true, false, bQiaA, bQijA, 1.0, 0.0);
    SharedTensor2d J_JAIK = std::make_shared<Tensor2d>("J (JA|IK)", naoccA, navirA, naoccA, naoccA);
    J_JAIK->sort(3214, J_IAJK, -1.0, 0.0);
    J_IAJK->axpy(J_JAIK, 1.0);
    J_JAIK.reset();
    SharedTensor2d G_IJAK = std::make_shared<Tensor2d>("G <IJ||AK>", naoccA, naoccA, navirA, naoccA);
    G_IJAK->sort(1324, J_IAJK, 1.0, 0.0);
    J_IAJK.reset();

    // form <IJ||AB>
    SharedTensor2d J_IAJB = std::make_shared<Tensor2d>("J (IA|JB)", naoccA, navirA, naoccA, navirA);
    J_IAJB->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    SharedTensor2d G_IJAB = std::make_shared<Tensor2d>("G <IJ||AB>", naoccA, naoccA, navirA, navirA);
    G_IJAB->sort(1324, J_IAJB, 1.0, 0.0);
    G_IJAB->sort(1342, J_IAJB, -1.0, 1.0);
    G_IJAB->write(psio_, PSIF_DFOCC_IJAB_AAAA);
    J_IAJB.reset();

    // malloc W[IJK](ABC)
    SharedTensor2d X_AAA = std::make_shared<Tensor2d>("X_AAA[IJK](A,CB)", navirA, navirA, navirA);
    SharedTensor2d Y_AAA = std::make_shared<Tensor2d>("Y_AAA[IJK](C,AB)", navirA, navirA, navirA);

    // malloc Wl[IJK](ABC)
    SharedTensor2d Xl_AAA = std::make_shared<Tensor2d>("Xl_AAA[IJK](A,CB)", navirA, navirA, navirA);
    SharedTensor2d Yl_AAA = std::make_shared<Tensor2d>("Yl_AAA[IJK](C,AB)", navirA, navirA, navirA);

    // malloc G[I](A,BC)
    G_I_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);
    SharedTensor2d G_J_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);
    SharedTensor2d G_K_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);

    // read T2AA
    SharedTensor2d t2AA = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    t2AA->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // read L2AA
    SharedTensor2d l2AA = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    l2AA->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // progress counter
    std::time_t stop, start = std::time(nullptr);
    long int ind = 0;
    double step_print = 10.0;
    double next_print = step_print;

    long int Nijk = naoccA * (naoccA - 1) * (naoccA - 2) / 6;
    outfile->Printf("\tNumber of ijk combinations for AAA: %i \n", Nijk);

    // N2 : main loop of AAA
    E_at = 0.0;
    double sumAAA = 0.0;
    for (long int i = 0; i < naoccA; i++) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);

        // read G[I](A,BC)
        G_I_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(i * navirA * navir2AA) * sizeof(double));
        for(long int j = 0; j < i; j++) {
            long int ij = ij_idxAA->get(i, j);
            double Dij = Di + FockA->get(j + nfrzc, j + nfrzc);

            // read G[J](A,BC)
            G_J_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(j * navirA * navir2AA) * sizeof(double));
            for (long int k = 0; k < j; k++) {
                long int kj = ij_idxAA->get(k, j);
                long int ik = ij_idxAA->get(i, k);

                // read G[K](A,BC)
                G_K_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(k * navirA * navir2AA) * sizeof(double));

                // X_AAA[IJK](A,CB) = \sum(E) t_IJ^AE <KE||CB>          (1)
                // X_AAA[IJK](A,CB) = \sum(E) T[IJ](A,E) G[K](E,CB)
                X_AAA->contract(false, false, navirA, navir2AA, navirA, t2AA, G_K_ABC,
                        (i * naoccA * navir2AA) + (j * navir2AA), 0, 1.0, 0.0);

                // X_AAA[IJK](A,CB) -= \sum(E) t_KJ^AE <IE||CB>         (2)
                // X_AAA[IJK](A,CB) -= \sum(E) T[KJ](A,E) G[I](E,CB)
                X_AAA->contract(false, false, navirA, navir2AA, navirA, t2AA, G_I_ABC,
                        (k * naoccA * navir2AA) + (j * navir2AA), 0, -1.0, 1.0);

                // X_AAA[IJK](A,CB) -= \sum(E) t_IK^AE <JE||CB>         (3)
                // X_AAA[IJK](A,CB) -= \sum(E) T[IK](A,E) G[J](E,CB)
                X_AAA->contract(false, false, navirA, navir2AA, navirA, t2AA, G_J_ABC,
                        (i * naoccA * navir2AA) + (k * navir2AA), 0, -1.0, 1.0);

                // Y_AAA[IJK](C,AB) = - \sum(M) t_IM^AB <KJ||CM         (4)
                // Y_AAA[IJK](C,AB) = - \sum(M) G[KJ](C,M) T[I](M,AB)
                Y_AAA->contract(false, false, navirA, navir2AA, naoccA, G_IJAK, t2AA,
                        (k * naoccA * navirA * naoccA) + (j * navirA * naoccA), (i * naoccA * navir2AA), -1.0, 0.0);

                // Y_AAA[IJK](C,AB) += \sum(M) t_JM^AB <KI||CM>         (5)
                // Y_AAA[IJK](C,AB) += \sum(M) G[KI](C,M) T[J](M,AB)
                Y_AAA->contract(false, false, navirA, navir2AA, naoccA, G_IJAK, t2AA,
                        (k * naoccA * navirA * naoccA) + (i * navirA * naoccA), (j * naoccA * navir2AA), 1.0, 1.0);

                // Y_AAA[IJK](C,AB) += \sum(M) t_KM^AB <IJ||CM>         (6)
                // Y_AAA[IJK](C,AB) += \sum(M) G[IJ](C,M) T[K](M,AB)
                Y_AAA->contract(false, false, navirA, navir2AA, naoccA, G_IJAK, t2AA,
                        (i * naoccA * navirA * naoccA) + (j * navirA * naoccA), (k * naoccA * navir2AA), 1.0, 1.0);

                // ==========================
                // === Asymmetric Triples ===
                // ==========================

                // Xl_AAA[IJK](A,CB) = \sum(E) l_IJ^AE <KE||CB>          (1)
                // Xl_AAA[IJK](A,CB) = \sum(E) L[IJ](A,E) G[K](E,CB)
                Xl_AAA->contract(false, false, navirA, navir2AA, navirA, l2AA, G_K_ABC,
                        (i * naoccA * navir2AA) + (j * navir2AA), 0, 1.0, 0.0);

                // Xl_AAA[IJK](A,CB) -= \sum(E) l_KJ^AE <IE||CB>         (2)
                // Xl_AAA[IJK](A,CB) -= \sum(E) L[KJ](A,E) G[I](E,CB)
                Xl_AAA->contract(false, false, navirA, navir2AA, navirA, l2AA, G_I_ABC,
                        (k * naoccA * navir2AA) + (j * navir2AA), 0, -1.0, 1.0);

                // Xl_AAA[IJK](A,CB) -= \sum(E) l_IK^AE <JE||CB>         (3)
                // Xl_AAA[IJK](A,CB) -= \sum(E) L[IK](A,E) G[J](E,CB)
                Xl_AAA->contract(false, false, navirA, navir2AA, navirA, l2AA, G_J_ABC,
                        (i * naoccA * navir2AA) + (k * navir2AA), 0, -1.0, 1.0);

                // Yl_AAA[IJK](C,AB) = - \sum(M) l_IM^AB <KJ||CM         (4)
                // Yl_AAA[IJK](C,AB) = - \sum(M) G[KJ](C,M) L[I](M,AB)
                Yl_AAA->contract(false, false, navirA, navir2AA, naoccA, G_IJAK, l2AA,
                        (k * naoccA * navirA * naoccA) + (j * navirA * naoccA), (i * naoccA * navir2AA), -1.0, 0.0);

                // Yl_AAA[IJK](C,AB) += \sum(M) l_JM^AB <KI||CM>         (5)
                // Yl_AAA[IJK](C,AB) += \sum(M) G[KI](C,M) L[J](M,AB)
                Yl_AAA->contract(false, false, navirA, navir2AA, naoccA, G_IJAK, l2AA,
                        (k * naoccA * navirA * naoccA) + (i * navirA * naoccA), (j * naoccA * navir2AA), 1.0, 1.0);

                // Yl_AAA[IJK](C,AB) += \sum(M) l_KM^AB <IJ||CM>         (6)
                // Yl_AAA[IJK](C,AB) += \sum(M) G[IJ](C,M) L[K](M,AB)
                Yl_AAA->contract(false, false, navirA, navir2AA, naoccA, G_IJAK, l2AA,
                        (i * naoccA * navirA * naoccA) + (j * navirA * naoccA), (k * naoccA * navir2AA), 1.0, 1.0);

                double Dijk = Dij + FockA->get(k + nfrzc, k + nfrzc);

                double Wijkabc, Wl_ijkabc, Vijkabc;
#pragma omp parallel for private(Wijkabc, Wl_ijkabc, Vijkabc) reduction(+ : sumAAA)
                for (long int a = 0; a < navirA; a++) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b < a; b++) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < b; c++) {
                            long int ac = ab_idxAA->get(a, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);

                            // W[IJK](ABC) = P(A/BC) X_AAA[IJK](ABC) + P(AB/C) Y_AAA[IJK](ABC)
                            Wijkabc = X_AAA->get(a, cb) - X_AAA->get(b, ca) - X_AAA->get(c, ab)
                                    + Y_AAA->get(c, ab) - Y_AAA->get(a, cb) - Y_AAA->get(b, ac);

                            // Wl[IJK](ABC) = P(A/BC) Xl_AAA[IJK](ABC) + P(AB/C) Yl_AAA[IJK](ABC)
                            Wl_ijkabc = Xl_AAA->get(a, cb) - Xl_AAA->get(b, ca) - Xl_AAA->get(c, ab)
                                      + Yl_AAA->get(c, ab) - Yl_AAA->get(a, cb) - Yl_AAA->get(b, ac);

                            // V_AAA = t_K^C <IJ||AB> + t_IJ^AB f_KC
                            //       - t_I^C <KJ||AB> - t_KJ^AB f_IC
                            //       - t_J^C <IK||AB> - t_IK^AB f_JC
                            //       - t_K^A <IJ||CB> + t_IJ^CB f_KA
                            //       + t_I^A <KJ||CB> - t_KJ^CB f_IA
                            //       + t_J^A <IK||CB> - t_IK^CB f_JA
                            //       - t_K^B <IJ||AC> + t_IJ^AC f_KB
                            //       + t_I^B <KJ||AC> - t_KJ^AC f_IB
                            //       + t_J^B <IK||AC> - t_IK^AC f_JB
                            //Vijkabc = t1A->get(k, c) * G_IJAB->get(ij, ab) + t2AA->get(ij, ab) * FockA->get(nfrzc + k, noccA + c)
                            //        - t1A->get(i, c) * G_IJAB->get(kj, ab) - t2AA->get(kj, ab) * FockA->get(nfrzc + i, noccA + c)
                            //        - t1A->get(j, c) * G_IJAB->get(ik, ab) - t2AA->get(ik, ab) * FockA->get(nfrzc + j, noccA + c)
                            //        - t1A->get(k, a) * G_IJAB->get(ij, cb) - t2AA->get(ij, cb) * FockA->get(nfrzc + k, noccA + a)
                            //        + t1A->get(i, a) * G_IJAB->get(kj, cb) + t2AA->get(kj, cb) * FockA->get(nfrzc + i, noccA + a)
                            //        + t1A->get(j, a) * G_IJAB->get(ik, cb) + t2AA->get(ik, cb) * FockA->get(nfrzc + j, noccA + a)
                            //        - t1A->get(k, b) * G_IJAB->get(ij, ac) - t2AA->get(ij, ac) * FockA->get(nfrzc + k, noccA + b)
                            //        + t1A->get(i, b) * G_IJAB->get(kj, ac) + t2AA->get(kj, ac) * FockA->get(nfrzc + i, noccA + b)
                            //        + t1A->get(j, b) * G_IJAB->get(ik, ac) + t2AA->get(ik, ac) * FockA->get(nfrzc + j, noccA + b);

                            // V_AAA = l_K^C <IJ||AB> + l_IJ^AB f_KC
                            //       - l_I^C <KJ||AB> - l_KJ^AB f_IC
                            //       - l_J^C <IK||AB> - l_IK^AB f_JC
                            //       - l_K^A <IJ||CB> + l_IJ^CB f_KA
                            //       + l_I^A <KJ||CB> - l_KJ^CB f_IA
                            //       + l_J^A <IK||CB> - l_IK^CB f_JA
                            //       - l_K^B <IJ||AC> + l_IJ^AC f_KB
                            //       + l_I^B <KJ||AC> - l_KJ^AC f_IB
                            //       + l_J^B <IK||AC> - l_IK^AC f_JB
                            Vijkabc = l1A->get(k, c) * G_IJAB->get(ij, ab) + l2AA->get(ij, ab) * FockA->get(nfrzc + k, noccA + c)
                                    - l1A->get(i, c) * G_IJAB->get(kj, ab) - l2AA->get(kj, ab) * FockA->get(nfrzc + i, noccA + c)
                                    - l1A->get(j, c) * G_IJAB->get(ik, ab) - l2AA->get(ik, ab) * FockA->get(nfrzc + j, noccA + c)
                                    - l1A->get(k, a) * G_IJAB->get(ij, cb) - l2AA->get(ij, cb) * FockA->get(nfrzc + k, noccA + a)
                                    + l1A->get(i, a) * G_IJAB->get(kj, cb) + l2AA->get(kj, cb) * FockA->get(nfrzc + i, noccA + a)
                                    + l1A->get(j, a) * G_IJAB->get(ik, cb) + l2AA->get(ik, cb) * FockA->get(nfrzc + j, noccA + a)
                                    - l1A->get(k, b) * G_IJAB->get(ij, ac) - l2AA->get(ij, ac) * FockA->get(nfrzc + k, noccA + b)
                                    + l1A->get(i, b) * G_IJAB->get(kj, ac) + l2AA->get(kj, ac) * FockA->get(nfrzc + i, noccA + b)
                                    + l1A->get(j, b) * G_IJAB->get(ik, ac) + l2AA->get(ik, ac) * FockA->get(nfrzc + j, noccA + b);

                            double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);

                            // N3 : energy of AAA
                            //sumAAA += (Wijkabc + Vijkabc) * Wl_ijkabc / Dijkabc;
                            sumAAA += (Wl_ijkabc + Vijkabc) * Wijkabc / Dijkabc;

                        } // c
                    } // b
                } // a

                // progress counter
                ind += 1;
                double percent = static_cast<double>(ind) / static_cast<double>(Nijk) * 100.0;
                if (percent >= next_print) {
                    stop = std::time(nullptr);
                    next_print += step_print;
                    outfile->Printf("              %5.1lf  %8d s\n", percent,
                                    static_cast<int>(stop) - static_cast<int>(start));
                }
            } // k
        } // j
    } // i

    outfile->Printf("\tAAA (AT) energy                    : % 20.14f\n", sumAAA);
    E_at += sumAAA;

    // reset all AAA things
    t2AA.reset();
    l2AA.reset();
    G_IJAB.reset();
    X_AAA.reset();
    Y_AAA.reset();
    G_I_ABC.reset();
    G_J_ABC.reset();
    G_K_ABC.reset();
    timer_off("CCSD(AT)-HM-AAA");

    timer_on("CCSD(AT)-HM-BBB");
    //==================================================
    //======================= BBB ======================
    //==================================================
    // N4 : beginning of BBB

    // form <ia||bc>
    // 'c' letter mean compact
    SharedTensor2d Jc_i_abc = std::make_shared<Tensor2d>("J[i] (a|b>=c)", navirB, ntri_abBB);
    //SharedTensor2d J_i_abc = std::make_shared<Tensor2d>("J[i] (a|bc)", navirB, navirB, navirB);
    SharedTensor2d G_i_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);
    SharedTensor2d bQabB_c = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ab)", nQ, ntri_abBB);
    bQabB_c->read(psio_, PSIF_DFOCC_INTS);
    SharedTensor2d biaQB = std::make_shared<Tensor2d>("B (ia|Q)", naoccB * navirB, nQ);
    biaQB->trans(bQiaB);

    for (long int i = 0; i < naoccB; i++) {
        bool flag = (i == 0) ? false : true;
        // [i](a|b>=c) = [i](a|Q) * (Q|a>=b)
        Jc_i_abc->contract(false, false, navirB, ntri_abBB, nQ, biaQB, bQabB_c, i * navirB * nQ, 0, 1.0, 0.0);

        // expand [i](a|b>=c) TO [i](a|bc)
        // [i]<b||ac> = [i]<b|ac> - [i]<b|ca>
        //            = [i](a|bc) - [i](c|ba)
        for (long int a = 0; a < navirB; a++) {
            for (long int b = 0; b < navirB; b++) {
                long int ab = ab_idxBB->get(a, b);
                long int ba = ab_idxBB->get(b, a);
                for (int c = 0; c < navirB; c++) {
                    long int ac = ab_idxBB->get(a, c);
                    long int ca = ab_idxBB->get(c, a);
                    long int bc = ab_idxBB->get(b, c);
                    //double val = Jc_i_abc->get(a, index2(b, c)); // for only expand
                    //J_i_abc->set(a, bc, val); // for only expand
                    //J_i_abc->set(b, ac, val); // for only sort
                    double val = Jc_i_abc->get(b, index2(a, c)) - Jc_i_abc->get(c, index2(b, a)); // for expand, sort and subtract (anti-symmetrization)
                    G_i_abc->set(a, bc, val); // for expand, sort and subtract (anti-symmetrization)
                } // c
            } // b
        } // a

        // write [i]<a||bc> to disk
        G_i_abc->mywrite(psio_, PSIF_DFOCC_IABC_BBBB, flag);
    } // i

    G_i_abc.reset();
    //J_i_abc.reset();
    Jc_i_abc.reset();

    // form <ij||ak>
    // (ia|jk) = bQiaB * bQijB
    // <ij||ak> = <ak||ij> = <ak|ij> - <ak|ji>
    //                     = (ai|kj) - (aj|ki)
    //                     = (ia|jk) - (ja|ik)
    // (ja|ik) = sort 3214 (ia|jk)
    SharedTensor2d J_iajk = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ia|jk)", naoccB, navirB, naoccB, naoccB);
    J_iajk->gemm(true, false, bQiaB, bQijB, 1.0, 0.0);
    SharedTensor2d J_jaik = std::make_shared<Tensor2d>("J (ja|ik)", naoccB, navirB, naoccB, naoccB);
    J_jaik->sort(3214, J_iajk, -1.0, 0.0);
    J_iajk->axpy(J_jaik, 1.0);
    J_jaik.reset();
    SharedTensor2d G_ijak = std::make_shared<Tensor2d>("G <ij||ak>", naoccB, naoccB, navirB, naoccB);
    G_ijak->sort(1324, J_iajk, 1.0, 0.0);
    J_iajk.reset();

    // form <ij||ab>
    SharedTensor2d J_iajb = std::make_shared<Tensor2d>("J (ia|jb)", naoccB, navirB, naoccB, navirB);
    J_iajb->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    SharedTensor2d G_ijab = std::make_shared<Tensor2d>("G <ij||ab>", naoccB, naoccB, navirB, navirB);
    G_ijab->sort(1324, J_iajb, 1.0, 0.0);
    G_ijab->sort(1342, J_iajb, -1.0, 1.0);
    G_ijab->write(psio_, PSIF_DFOCC_IJAB_BBBB);
    J_iajb.reset();

    // malloc W[ijk](abc)
    SharedTensor2d X_BBB = std::make_shared<Tensor2d>("X[ijk](a,cb)", navirB, navirB, navirB);
    SharedTensor2d Y_BBB = std::make_shared<Tensor2d>("Y[ijk](c,ab)", navirB, navirB, navirB);

    // malloc Wl[ijk](abc)
    SharedTensor2d Xl_BBB = std::make_shared<Tensor2d>("Xl[ijk](a,cb)", navirB, navirB, navirB);
    SharedTensor2d Yl_BBB = std::make_shared<Tensor2d>("Yl[ijk](c,ab)", navirB, navirB, navirB);

    // read T2BB
    SharedTensor2d t2BB = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    t2BB->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // read L2BB
    SharedTensor2d l2BB = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    l2BB->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // malloc G[i](a,bc)
    G_i_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);
    SharedTensor2d G_j_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);
    SharedTensor2d G_k_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);

    // progress counter
    start = std::time(nullptr);
    stop = std::time(nullptr);
    ind = 0;
    next_print = step_print;

    Nijk = naoccB * (naoccB - 1) * (naoccB - 2) / 6;
    outfile->Printf("\n\tNumber of ijk combinations for BBB: %i \n", Nijk);

    // N5 : main loop of BBB
    double sumBBB = 0.0;
    for (long int i = 0; i < naoccB; i++) {
        double Di = FockB->get(i + nfrzc, i + nfrzc);

        // read G[i](a,bc)
        G_i_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(i * navirB * navir2BB) * sizeof(double));
        for(long int j = 0; j < i; j++) {
            long int ij = ij_idxBB->get(i, j);
            double Dij = Di + FockB->get(j + nfrzc, j + nfrzc);

            // read G[j](a,bc)
            G_j_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(j * navirB * navir2BB) * sizeof(double));
            for (long int k = 0; k < j; k++) {
                long int kj = ij_idxBB->get(k, j);
                long int ik = ij_idxBB->get(i, k);

                // read G[k](a,bc)
                G_k_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(k * navirB * navir2BB) * sizeof(double));

                // X_BBB[ijk](a,cb) = \sum(e) t_ij^ae <ke||cb>          (1)
                // X_BBB[ijk](a,cb) = \sum(e) T[ij](a,e) G[k](e,cb)
                X_BBB->contract(false, false, navirB, navir2BB, navirB, t2BB, G_k_abc,
                        (i * naoccB * navir2BB) + (j * navir2BB), 0, 1.0, 0.0);

                // X_BBB[ijk](a,cb) -= \sum(e) t_kj^ae <ie||cb>         (2)
                // X_BBB[ijk](a,cb) -= \sum(e) T[kj](a,e) G[i](e,cb)
                X_BBB->contract(false, false, navirB, navir2BB, navirB, t2BB, G_i_abc,
                        (k * naoccB * navir2BB) + (j * navir2BB), 0, -1.0, 1.0);

                // X_BBB[ijk](a,cb) -= \sum(e) t_ik^ae <je||cb>         (3)
                // X_BBB[ijk](a,cb) -= \sum(e) T[ik](a,e) G[j](e,cb)
                X_BBB->contract(false, false, navirB, navir2BB, navirB, t2BB, G_j_abc,
                        (i * naoccB * navir2BB) + (k * navir2BB), 0, -1.0, 1.0);

                // Y_BBB[ijk](c,ab) = - \sum(m) t_im^ab <kj||cm>        (4)
                // Y_BBB[ijk](c,ab) = - \sum(m) G[kj](c,m) T[i](m,ab)
                Y_BBB->contract(false, false, navirB, navir2BB, naoccB, G_ijak, t2BB,
                        (k * naoccB * navirB * naoccB) + (j * navirB * naoccB), (i * naoccB * navir2BB), -1.0, 0.0);

                // Y_BBB[ijk](c,ab) += \sum(m) t_jm^ab <ki||cm>         (5)
                // Y_BBB[ijk](c,ab) += \sum(m) G[ki](c,m) T[j](m,ab)
                Y_BBB->contract(false, false, navirB, navir2BB, naoccB, G_ijak, t2BB,
                        (k * naoccB * navirB * naoccB) + (i * navirB * naoccB), (j * naoccB * navir2BB), 1.0, 1.0);

                // Y_BBB[ijk](c,ab) += \sum(m) t_km^ab <ij||cm>         (6)
                // Y_BBB[ijk](c,ab) += \sum(m) G[ij](c,m) T[k](m,ab)
                Y_BBB->contract(false, false, navirB, navir2BB, naoccB, G_ijak, t2BB,
                        (i * naoccB * navirB * naoccB) + (j * navirB * naoccB), (k * naoccB * navir2BB), 1.0, 1.0);

                // ==========================
                // === Asymmetric Triples ===
                // ==========================

                // Xl_BBB[ijk](a,cb) = \sum(e) l_ij^ae <ke||cb>          (1)
                // Xl_BBB[ijk](a,cb) = \sum(e) L[ij](a,e) G[k](e,cb)
                Xl_BBB->contract(false, false, navirB, navir2BB, navirB, l2BB, G_k_abc,
                        (i * naoccB * navir2BB) + (j * navir2BB), 0, 1.0, 0.0);

                // Xl_BBB[ijk](a,cb) -= \sum(e) l_kj^ae <ie||cb>         (2)
                // Xl_BBB[ijk](a,cb) -= \sum(e) L[kj](a,e) G[i](e,cb)
                Xl_BBB->contract(false, false, navirB, navir2BB, navirB, l2BB, G_i_abc,
                        (k * naoccB * navir2BB) + (j * navir2BB), 0, -1.0, 1.0);

                // Xl_BBB[ijk](a,cb) -= \sum(e) l_ik^ae <je||cb>         (3)
                // Xl_BBB[ijk](a,cb) -= \sum(e) L[ik](a,e) G[j](e,cb)
                Xl_BBB->contract(false, false, navirB, navir2BB, navirB, l2BB, G_j_abc,
                        (i * naoccB * navir2BB) + (k * navir2BB), 0, -1.0, 1.0);

                // Yl_BBB[ijk](c,ab) = - \sum(m) l_im^ab <kj||cm>        (4)
                // Yl_BBB[ijk](c,ab) = - \sum(m) G[kj](c,m) L[i](m,ab)
                Yl_BBB->contract(false, false, navirB, navir2BB, naoccB, G_ijak, l2BB,
                        (k * naoccB * navirB * naoccB) + (j * navirB * naoccB), (i * naoccB * navir2BB), -1.0, 0.0);

                // Yl_BBB[ijk](c,ab) += \sum(m) l_jm^ab <ki||cm>         (5)
                // Yl_BBB[ijk](c,ab) += \sum(m) G[ki](c,m) T[j](m,ab)
                Yl_BBB->contract(false, false, navirB, navir2BB, naoccB, G_ijak, l2BB,
                        (k * naoccB * navirB * naoccB) + (i * navirB * naoccB), (j * naoccB * navir2BB), 1.0, 1.0);

                // Yl_BBB[ijk](c,ab) += \sum(m) l_km^ab <ij||cm>         (6)
                // Yl_BBB[ijk](c,ab) += \sum(m) G[ij](c,m) L[k](m,ab)
                Yl_BBB->contract(false, false, navirB, navir2BB, naoccB, G_ijak, l2BB,
                        (i * naoccB * navirB * naoccB) + (j * navirB * naoccB), (k * naoccB * navir2BB), 1.0, 1.0);

                double Dijk = Dij + FockB->get(k + nfrzc, k + nfrzc);

                double Wijkabc, Wl_ijkabc, Vijkabc;
#pragma omp parallel for private(Wijkabc, Wl_ijkabc ,Vijkabc) reduction(+ : sumBBB)
                for (long int a = 0; a < navirB; a++) {
                    double Dijka = Dijk - FockB->get(a + noccB, a + noccB);
                    for (long int b = 0; b < a; b++) {
                        double Dijkab = Dijka - FockB->get(b + noccB, b + noccB);
                        long int ab = ab_idxBB->get(a, b);
                        for (long int c = 0; c < b; c++) {
                            long int ac = ab_idxBB->get(a, c);
                            long int ca = ab_idxBB->get(c, a);
                            long int cb = ab_idxBB->get(c, b);

                            // W[ijk](abc) = P(a/bc) X_BBB[ijk](abc) + P(ab/c) Y_BBB[ijk](abc)
                            Wijkabc = X_BBB->get(a, cb) - X_BBB->get(b, ca) - X_BBB->get(c, ab)
                                    + Y_BBB->get(c, ab) - Y_BBB->get(a, cb) - Y_BBB->get(b, ac);

                            // Wl[ijk](abc) = P(a/bc) Xl_BBB[ijk](abc) + P(ab/c) Yl_BBB[ijk](abc)
                            Wl_ijkabc = Xl_BBB->get(a, cb) - Xl_BBB->get(b, ca) - Xl_BBB->get(c, ab)
                                      + Yl_BBB->get(c, ab) - Yl_BBB->get(a, cb) - Yl_BBB->get(b, ac);

                            // V_BBB = t_k^c <ij||ab> + t_ij^ab f_kc
                            //       - t_i^c <kj||ab> - t_kj^ab f_ic
                            //       - t_j^c <ik||ab> - t_ik^ab f_jc
                            //       - t_k^a <ij||cb> + t_ij^cb f_ka
                            //       + t_i^a <kj||cb> - t_kj^cb f_ia
                            //       + t_j^a <ik||cb> - t_ik^cb f_ja
                            //       - t_k^b <ij||ac> + t_ij^ac f_kb
                            //       + t_i^b <kj||ac> - t_kj^ac f_ib
                            //       + t_j^b <ik||ac> - t_ik^ac f_jb
                            //Vijkabc = t1B->get(k, c) * G_ijab->get(ij, ab) + t2BB->get(ij, ab) * FockB->get(nfrzc + k, noccB + c)
                            //        - t1B->get(i, c) * G_ijab->get(kj, ab) - t2BB->get(kj, ab) * FockB->get(nfrzc + i, noccB + c)
                            //        - t1B->get(j, c) * G_ijab->get(ik, ab) - t2BB->get(ik, ab) * FockB->get(nfrzc + j, noccB + c)
                            //        - t1B->get(k, a) * G_ijab->get(ij, cb) - t2BB->get(ij, cb) * FockB->get(nfrzc + k, noccB + a)
                            //        + t1B->get(i, a) * G_ijab->get(kj, cb) + t2BB->get(kj, cb) * FockB->get(nfrzc + i, noccB + a)
                            //        + t1B->get(j, a) * G_ijab->get(ik, cb) + t2BB->get(ik, cb) * FockB->get(nfrzc + j, noccB + a)
                            //        - t1B->get(k, b) * G_ijab->get(ij, ac) - t2BB->get(ij, ac) * FockB->get(nfrzc + k, noccB + b)
                            //        + t1B->get(i, b) * G_ijab->get(kj, ac) + t2BB->get(kj, ac) * FockB->get(nfrzc + i, noccB + b)
                            //        + t1B->get(j, b) * G_ijab->get(ik, ac) + t2BB->get(ik, ac) * FockB->get(nfrzc + j, noccB + b);

                            // V_BBB = l_k^c <ij||ab> + l_ij^ab f_kc
                            //       - l_i^c <kj||ab> - l_kj^ab f_ic
                            //       - l_j^c <ik||ab> - l_ik^ab f_jc
                            //       - l_k^a <ij||cb> + l_ij^cb f_ka
                            //       + l_i^a <kj||cb> - l_kj^cb f_ia
                            //       + l_j^a <ik||cb> - l_ik^cb f_ja
                            //       - l_k^b <ij||ac> + l_ij^ac f_kb
                            //       + l_i^b <kj||ac> - l_kj^ac f_ib
                            //       + l_j^b <ik||ac> - l_ik^ac f_jb
                            Vijkabc = l1B->get(k, c) * G_ijab->get(ij, ab) + l2BB->get(ij, ab) * FockB->get(nfrzc + k, noccB + c)
                                    - l1B->get(i, c) * G_ijab->get(kj, ab) - l2BB->get(kj, ab) * FockB->get(nfrzc + i, noccB + c)
                                    - l1B->get(j, c) * G_ijab->get(ik, ab) - l2BB->get(ik, ab) * FockB->get(nfrzc + j, noccB + c)
                                    - l1B->get(k, a) * G_ijab->get(ij, cb) - l2BB->get(ij, cb) * FockB->get(nfrzc + k, noccB + a)
                                    + l1B->get(i, a) * G_ijab->get(kj, cb) + l2BB->get(kj, cb) * FockB->get(nfrzc + i, noccB + a)
                                    + l1B->get(j, a) * G_ijab->get(ik, cb) + l2BB->get(ik, cb) * FockB->get(nfrzc + j, noccB + a)
                                    - l1B->get(k, b) * G_ijab->get(ij, ac) - l2BB->get(ij, ac) * FockB->get(nfrzc + k, noccB + b)
                                    + l1B->get(i, b) * G_ijab->get(kj, ac) + l2BB->get(kj, ac) * FockB->get(nfrzc + i, noccB + b)
                                    + l1B->get(j, b) * G_ijab->get(ik, ac) + l2BB->get(ik, ac) * FockB->get(nfrzc + j, noccB + b);

                            double Dijkabc = Dijkab - FockB->get(c + noccB, c + noccB);

                            // N6 : energy of BBB
                            //sumBBB += (Wijkabc + Vijkabc) * Wl_ijkabc / Dijkabc;
                            sumBBB += (Wl_ijkabc + Vijkabc) * Wijkabc / Dijkabc;

                        } // c
                    } // b
                } // a

                // progress counter
                ind += 1;
                double percent = static_cast<double>(ind) / static_cast<double>(Nijk) * 100.0;
                if (percent >= next_print) {
                    stop = std::time(nullptr);
                    next_print += step_print;
                    outfile->Printf("              %5.1lf  %8d s\n", percent,
                                    static_cast<int>(stop) - static_cast<int>(start));
                }
            } // k
        } // j
    } // i

    outfile->Printf("\tBBB (AT) energy                    : % 20.14f\n", sumBBB);
    E_at += sumBBB;

    // reset all BBB things
    t2BB.reset();
    l2BB.reset();
    G_ijab.reset();
    X_BBB.reset();
    Y_BBB.reset();
    Xl_BBB.reset();
    Yl_BBB.reset();
    G_i_abc.reset();
    timer_off("CCSD(AT)-HM-BBB");

    //==================================================
    //======================= AAB ======================
    //==================================================
    // N7 : beginning of AAB
    timer_on("CCSD(AT)-HM-AAB");

    // form <iA|bC>
    SharedTensor2d Jc_i_aBC = std::make_shared<Tensor2d>("J[i] (a|B>=C)", navirB, ntri_abAA);
    SharedTensor2d K_i_AbC = std::make_shared<Tensor2d>("K[i] <A|bC>", navirA, navirB, navirA);
    bQabA_c = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    bQabA_c->read(psio_, PSIF_DFOCC_INTS);
    biaQB->trans(bQiaB);
    for (long int i = 0; i < naoccB; i++) {
        bool flag = (i == 0) ? false : true;
        // [i](a|B>=C) = [i](a|Q) * (Q|A>=B)
        Jc_i_aBC->contract(false, false, navirB, ntri_abAA, nQ, biaQB, bQabA_c, i * navirB * nQ, 0, 1.0, 0.0);

        // expand [i](a|B>=C) to [i](a|BC)
        // [i]<A|bC> = sort 213 [i](a|BC)
        for (long int a = 0; a < navirB; a++) {
            for (long int b = 0; b < navirA; b++) {
                for (long int c = 0; c < navirA; c++) {
                    long int ac = ab_idxBA->get(a, c);
                    long int bc = ab_idxAA->get(b, c);
                    double val = Jc_i_aBC->get(a, index2(b, c));
                    K_i_AbC->set(b, ac, val);
                } // c
            } // b
        } // a

        // write [i]<A|bC> to disk
        K_i_AbC->mywrite(psio_, PSIF_DFOCC_IABC_BABA, flag);
    } // i

    K_i_AbC.reset();
    Jc_i_aBC.reset();
    bQabA_c.reset();
    biaQB.reset();

    // form <Ia|Bc>
    SharedTensor2d K_I_aBc = std::make_shared<Tensor2d>("K[I] <a|Bc>", navirB, navirA, navirB);
    SharedTensor2d Jc_I_Abc = std::make_shared<Tensor2d>("J[I] (A|b>=c)", navirA, ntri_abBB);
    bQabB_c->read(psio_, PSIF_DFOCC_INTS);
    biaQA->trans(bQiaA);

    for (long int i = 0; i < naoccA; i++) {
        bool flag = (i == 0) ? false : true;
        // [I](A|b>=c) = [I](A|Q) * (Q|a>=b)
        Jc_I_Abc->contract(false, false, navirA, ntri_abBB, nQ, biaQA, bQabB_c, i * navirA * nQ, 0, 1.0, 0.0);

        // expand [I](A|b>=c) to [I](A|bc)
        // [I]<a|Bc> = sort 213 [I](A|bc)
        for (long int a = 0; a < navirA; a++) {
            for (long int b = 0; b < navirB; b++) {
                for (long int c = 0; c < navirB; c++) {
                    long int ac = ab_idxAB->get(a, c);
                    long int bc = ab_idxBB->get(b, c);
                    double val = Jc_I_Abc->get(a, index2(b, c));
                    K_I_aBc->set(b, ac, val);
                } // c
            } // b
        } // a

        // write [I]<a|Bc> to disk
        K_I_aBc->mywrite(psio_, PSIF_DFOCC_IABC_ABAB, flag);
    } // i

    K_I_aBc.reset();
    Jc_I_Abc.reset();
    bQabB_c.reset();
    biaQA.reset();

    // form <Ij|Ka>
    SharedTensor2d J_IJka = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|ka)", naoccA, naoccA, naoccB, navirB);
    J_IJka->gemm(true, false, bQijA, bQiaB, 1.0, 0.0);
    SharedTensor2d K_IjKa = std::make_shared<Tensor2d>("K <Ij|Ka>", naoccA, naoccB, naoccA, navirB);
    K_IjKa->sort(1324, J_IJka, 1.0, 0.0);
    J_IJka.reset();

    // form <Ij|Ak>
    SharedTensor2d J_IAjk = std::make_shared<Tensor2d>("J (IA|jk)", naoccA, navirA, naoccB, naoccB);
    J_IAjk->gemm(true, false, bQiaA, bQijB, 1.0, 0.0);
    SharedTensor2d K_IjAk = std::make_shared<Tensor2d>("K <Ij|Ak>", naoccA, naoccB, navirA, naoccB);
    K_IjAk->sort(1324, J_IAjk, 1.0, 0.0);
    J_IAjk.reset();

    // form <Ij|Ab>
    SharedTensor2d J_IAjb = std::make_shared<Tensor2d>("J (IA|jb)", naoccA, navirA, naoccB, navirB);
    J_IAjb->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    SharedTensor2d K_IjAb = std::make_shared<Tensor2d>("K <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    K_IjAb->sort(1324, J_IAjb, 1.0, 0.0);
    K_IjAb->write(psio_, PSIF_DFOCC_IJAB_ABAB);
    J_IAjb.reset();

    SharedTensor2d W_AAB;
    SharedTensor2d Wl_AAB;
    SharedTensor2d Waux;

    // read T2AA
    t2AA = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    t2AA->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // read T2AB
    SharedTensor2d t2AB = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    t2AB->read(psio_, PSIF_DFOCC_AMPS);

    // form T2BA
    SharedTensor2d t2BA = std::make_shared<Tensor2d>("T2 <iJ|aB>", naoccB, naoccA, navirB, navirA);
    t2BA->sort(2143, t2AB, 1.0, 0.0);

    // read L2AA
    l2AA = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    l2AA->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // read L2AB
    SharedTensor2d l2AB = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    l2AB->read(psio_, PSIF_DFOCC_AMPS);

    // form L2BA
    SharedTensor2d l2BA = std::make_shared<Tensor2d>("L2 <iJ|aB>", naoccB, naoccA, navirB, navirA);
    l2BA->sort(2143, l2AB, 1.0, 0.0);

    G_IJAB = std::make_shared<Tensor2d>("G <IJ||AB>", naoccA, naoccA, navirA, navirA);
    G_IJAB->read(psio_, PSIF_DFOCC_IJAB_AAAA);

    // progress counter
    start = std::time(nullptr);
    stop = std::time(nullptr);
    ind = 0;
    next_print = step_print;

    Nijk = naoccA * (naoccA - 1) * (naoccB) / 2;
    outfile->Printf("\n\tNumber of ijk combinations for AAB: %i \n", Nijk);

    // N8 : main loop of AAB
    double sumAAB = 0.0;
    for (long int i = 0; i < naoccA; i++) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);

        for(long int j = 0; j < i; j++) {
            long int ij = ij_idxAA->get(i, j);
            double Dij = Di + FockA->get(j + nfrzc, j + nfrzc);

            for (long int k = 0; k < naoccB; k++) {
                long int kj = ij_idxBA->get(k, j);
                long int ik = ij_idxAB->get(i, k);
                long int jk = ij_idxAB->get(j, k);

                // read K[k](A,bC)
                K_i_AbC = std::make_shared<Tensor2d>("K[i] <A|bC>", navirA, navirB, navirA);
                K_i_AbC->myread(psio_, PSIF_DFOCC_IABC_BABA, (size_t)(k * navirB * navir2AA) * sizeof(double));

                Waux = std::make_shared<Tensor2d>("X[IJk](A, cB)", navirA, navirB, navirA);
                // Waux[IJk](A,cB) = \sum(E) t_IJ^AE <kE|cB>          (1)
                // Waux[IJk](A,cB) = \sum(E) T[IJ](A,E) K[k](E,cB)
                Waux->contract(false, false, navirA, navirB * navirA, navirA, t2AA, K_i_AbC,
                        (i * naoccA * navir2AA) + (j * navir2AA), 0, 1.0, 0.0);
                K_i_AbC.reset();

                W_AAB = std::make_shared<Tensor2d>("X[IJk](AB, c)", navirA * navirA, navirB);

                // Waux(A,cB) sort 132 --> Waux(AB,c)
                // Waux(A,cB) sort 312 --> Waux(BA,c)
                // W_AAB(AB,c) += Waux(AB,c) - Waux(BA,c)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int B = 0; B < navirA; ++B) {
                        W_AAB->axpy((size_t)navirB, A * navirA * navirB + B, navirA, Waux, A * navirA * navirB + B * navirB, 1, 1.0);
                        W_AAB->axpy((size_t)navirB, A * navirA * navirB + B, navirA, Waux, B * navirA * navirB + A * navirB, 1, -1.0);
                    }
                }

                // Waux[IJk](B,cA) = \sum(M) t_kM^cA <IJ||BM>         (2)
                // Waux[IJk](B,cA) = \sum(M) G[IJ](B,M) T[k](M,cA)
                Waux->contract(false, false, navirA, navirB * navirA, naoccA, G_IJAK, t2BA,
                        (i * naoccA * navirA * naoccA) + (j * navirA * naoccA), (k * naoccA * navirB * navirA), 1.0, 0.0);

                // Waux(B,cA) sort 312 --> Waux(AB,c)
                // Waux(B,cA) sort 132 --> Waux(BA,c)
                // W_AAB(AB,c) += Waux(AB,c) - Waux(BA,c)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int B = 0; B < navirA; ++B) {
                        W_AAB->axpy((size_t)navirB, B * navirA * navirB + A, navirA, Waux, A * navirA * navirB + B * navirB, 1, 1.0);
                        W_AAB->axpy((size_t)navirB, B * navirA * navirB + A, navirA, Waux, B * navirA * navirB + A * navirB, 1, -1.0);
                    }
                }
                Waux.reset();

                // read G[I](A,BC)
                G_I_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);
                G_I_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(i * navirA * navir2AA) * sizeof(double));

                // W_AAB[IJk](AB,c) = \sum(E) t_Jk^Ec <IE||AB>         (3)
                // W_AAB[IJk](AB,c) = \sum(E) G[I](E,AB) T[Jk](E,c)
                W_AAB->contract(true, false, navir2AA, navirB, navirA, G_I_ABC, t2AB,
                        0, (j * naoccB * navirA * navirB) + (k * navirA * navirB), 1.0, 1.0);

                // read G[J](A,BC)
                G_I_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(j * navirA * navir2AA) * sizeof(double));

                // W_AAB[IJk](AB,c) -= \sum(E) t_Ik^Ec <JE||AB>        (4)
                // W_AAB[IJk](AB,c) -= \sum(E) G[J](E,AB) T[Ik](E,c)
                W_AAB->contract(true, false, navir2AA, navirB, navirA, G_I_ABC, t2AB,
                        0, (i * naoccB * navirA * navirB) + (k * navirA * navirB), -1.0, 1.0);
                G_I_ABC.reset();

                // W_AAB[IJk](AB,c) += \sum(M) t_JM^AB <Ik|Mc>         (5)
                // W_AAB[IJk](AB,c) += \sum(M) T[J](M,AB) K[Ik](M,c)
                W_AAB->contract(true, false, navir2AA, navirB, naoccA, t2AA, K_IjKa,
                        (j * naoccA * navir2AA), (i * naoccB * naoccA * navirB) + (k * naoccA * navirB), 1.0, 1.0);

                // W_AAB[IJk](AB,c) -= \sum(M) t_IM^AB <Jk|Mc>         (6)
                // W_AAB[IJk](AB,c) -= \sum(M) T[I](M,AB) K[Jk](M,c)
                W_AAB->contract(true, false, navir2AA, navirB, naoccA, t2AA, K_IjKa,
                        (i * naoccA * navir2AA), (j * naoccB * naoccA * navirB) + (k * naoccA * navirB), -1.0, 1.0);

                // read K[J](a,Bc)
                K_I_aBc = std::make_shared<Tensor2d>("K[I] <a|Bc>", navirB, navirA, navirB);
                K_I_aBc->myread(psio_, PSIF_DFOCC_IABC_ABAB, (size_t)(j * navirB * navirA * navirB) * sizeof(double));

                Waux = std::make_shared<Tensor2d>("X[IJk](A, Bc)", navirA, navirA, navirB);

                // Waux[IJk](A,Bc) = \sum(e) t_Ik^Ae <Je|Bc>          (7)
                // Waux[IJk](A,Bc) = \sum(e) T[Ik](A,e) K[J](e,Bc)
                Waux->contract(false, false, navirA, navirA * navirB, navirB, t2AB, K_I_aBc,
                        (i * naoccB * navirA * navirB) + (k * navirA * navirB), 0, 1.0, 0.0);

                // read K[I](a,Bc)
                K_I_aBc->myread(psio_, PSIF_DFOCC_IABC_ABAB, (size_t)(i * navirB * navirA * navirB) * sizeof(double));

                // Waux[IJk](A,Bc) -= \sum(e) t_Jk^Ae <Ie|Bc>         (8)
                // Waux[IJk](A,Bc) -= \sum(e) T[Jk](A,e) K[I](e,Bc)
                Waux->contract(false, false, navirA, navirA * navirB, navirB, t2AB, K_I_aBc,
                        (j * naoccB * navirA * navirB) + (k * navirA * navirB), 0, -1.0, 1.0);
                K_I_aBc.reset();

                // Waux[IJk](A,Bc) += \sum(m) t_Im^Bc <Jk|Am>         (9)
                // Waux[IJk](A,Bc) += \sum(m) K[Jk](A,m) T[I](m,Bc)
                Waux->contract(false, false, navirA, navirA * navirB, naoccB, K_IjAk, t2AB,
                        (j * naoccB * navirA * naoccB) + (k * navirA * naoccB), (i * naoccB * navirA * navirB), 1.0, 1.0);

                // Waux[IJk](A,Bc) -= \sum(m) t_Jm^Bc <Ik|Am>         (10)
                // Waux[IJk](A,Bc) -= \sum(m) K[Ik](A,m) T[J](m,Bc)
                Waux->contract(false, false, navirA, navirA * navirB, naoccB, K_IjAk, t2AB,
                        (i * naoccB * navirA * naoccB) + (k * navirA * naoccB), (j * naoccB * navirA * navirB), -1.0, 1.0);

                // W_AAB(AB,c) += Waux(AB,c)
                W_AAB->axpy(Waux, 1.0);
                // Waux(A,Bc) sort 213 --> Waux(BA,c)
                // W_AAB(AB,c) -= Waux(BA,c)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int B = 0; B < navirA; ++B) {
                        W_AAB->axpy((size_t)navirB, A * navirA * navirB + B * navirB, 1, Waux, B * navirA * navirB + A * navirB, 1, -1.0);
                    }
                }
                Waux.reset();

                // ==========================
                // === Asymmetric Triples ===
                // ==========================

                // read K[k](A,bC)
                K_i_AbC = std::make_shared<Tensor2d>("K[i] <A|bC>", navirA, navirB, navirA);
                K_i_AbC->myread(psio_, PSIF_DFOCC_IABC_BABA, (size_t)(k * navirB * navir2AA) * sizeof(double));

                Waux = std::make_shared<Tensor2d>("X[IJk](A, cB)", navirA, navirB, navirA);
                // Waux[IJk](A,cB) = \sum(E) l_IJ^AE <kE|cB>          (1)
                // Waux[IJk](A,cB) = \sum(E) L[IJ](A,E) K[k](E,cB)
                Waux->contract(false, false, navirA, navirB * navirA, navirA, l2AA, K_i_AbC,
                        (i * naoccA * navir2AA) + (j * navir2AA), 0, 1.0, 0.0);
                K_i_AbC.reset();

                Wl_AAB = std::make_shared<Tensor2d>("X[IJk](AB, c)", navirA * navirA, navirB);

                // Waux(A,cB) sort 132 --> Waux(AB,c)
                // Waux(A,cB) sort 312 --> Waux(BA,c)
                // Wl_AAB(AB,c) += Waux(AB,c) - Waux(BA,c)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int B = 0; B < navirA; ++B) {
                        Wl_AAB->axpy((size_t)navirB, A * navirA * navirB + B, navirA, Waux, A * navirA * navirB + B * navirB, 1, 1.0);
                        Wl_AAB->axpy((size_t)navirB, A * navirA * navirB + B, navirA, Waux, B * navirA * navirB + A * navirB, 1, -1.0);
                    }
                }

                // Waux[IJk](B,cA) = \sum(M) l_kM^cA <IJ||BM>         (2)
                // Waux[IJk](B,cA) = \sum(M) G[IJ](B,M) L[k](M,cA)
                Waux->contract(false, false, navirA, navirB * navirA, naoccA, G_IJAK, l2BA,
                        (i * naoccA * navirA * naoccA) + (j * navirA * naoccA), (k * naoccA * navirB * navirA), 1.0, 0.0);

                // Waux(B,cA) sort 312 --> Waux(AB,c)
                // Waux(B,cA) sort 132 --> Waux(BA,c)
                // Wl_AAB(AB,c) += Waux(AB,c) - Waux(BA,c)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int B = 0; B < navirA; ++B) {
                        Wl_AAB->axpy((size_t)navirB, B * navirA * navirB + A, navirA, Waux, A * navirA * navirB + B * navirB, 1, 1.0);
                        Wl_AAB->axpy((size_t)navirB, B * navirA * navirB + A, navirA, Waux, B * navirA * navirB + A * navirB, 1, -1.0);
                    }
                }
                Waux.reset();

                // read G[I](A,BC)
                G_I_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);
                G_I_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(i * navirA * navir2AA) * sizeof(double));

                // Wl_AAB[IJk](AB,c) = \sum(E) l_Jk^Ec <IE||AB>         (3)
                // Wl_AAB[IJk](AB,c) = \sum(E) G[I](E,AB) L[Jk](E,c)
                Wl_AAB->contract(true, false, navir2AA, navirB, navirA, G_I_ABC, l2AB,
                        0, (j * naoccB * navirA * navirB) + (k * navirA * navirB), 1.0, 1.0);

                // read G[J](A,BC)
                G_I_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(j * navirA * navir2AA) * sizeof(double));

                // Wl_AAB[IJk](AB,c) -= \sum(E) l_Ik^Ec <JE||AB>        (4)
                // Wl_AAB[IJk](AB,c) -= \sum(E) G[J](E,AB) L[Ik](E,c)
                Wl_AAB->contract(true, false, navir2AA, navirB, navirA, G_I_ABC, l2AB,
                        0, (i * naoccB * navirA * navirB) + (k * navirA * navirB), -1.0, 1.0);
                G_I_ABC.reset();

                // Wl_AAB[IJk](AB,c) += \sum(M) l_JM^AB <Ik|Mc>         (5)
                // Wl_AAB[IJk](AB,c) += \sum(M) L[J](M,AB) K[Ik](M,c)
                Wl_AAB->contract(true, false, navir2AA, navirB, naoccA, l2AA, K_IjKa,
                        (j * naoccA * navir2AA), (i * naoccB * naoccA * navirB) + (k * naoccA * navirB), 1.0, 1.0);

                // Wl_AAB[IJk](AB,c) -= \sum(M) l_IM^AB <Jk|Mc>         (6)
                // Wl_AAB[IJk](AB,c) -= \sum(M) L[I](M,AB) K[Jk](M,c)
                Wl_AAB->contract(true, false, navir2AA, navirB, naoccA, l2AA, K_IjKa,
                        (i * naoccA * navir2AA), (j * naoccB * naoccA * navirB) + (k * naoccA * navirB), -1.0, 1.0);

                // read K[J](a,Bc)
                K_I_aBc = std::make_shared<Tensor2d>("K[I] <a|Bc>", navirB, navirA, navirB);
                K_I_aBc->myread(psio_, PSIF_DFOCC_IABC_ABAB, (size_t)(j * navirB * navirA * navirB) * sizeof(double));

                Waux = std::make_shared<Tensor2d>("X[IJk](A, Bc)", navirA, navirA, navirB);

                // Waux[IJk](A,Bc) = \sum(e) l_Ik^Ae <Je|Bc>          (7)
                // Waux[IJk](A,Bc) = \sum(e) L[Ik](A,e) K[J](e,Bc)
                Waux->contract(false, false, navirA, navirA * navirB, navirB, l2AB, K_I_aBc,
                        (i * naoccB * navirA * navirB) + (k * navirA * navirB), 0, 1.0, 0.0);

                // read K[I](a,Bc)
                K_I_aBc->myread(psio_, PSIF_DFOCC_IABC_ABAB, (size_t)(i * navirB * navirA * navirB) * sizeof(double));

                // Waux[IJk](A,Bc) -= \sum(e) l_Jk^Ae <Ie|Bc>         (8)
                // Waux[IJk](A,Bc) -= \sum(e) L[Jk](A,e) K[I](e,Bc)
                Waux->contract(false, false, navirA, navirA * navirB, navirB, l2AB, K_I_aBc,
                        (j * naoccB * navirA * navirB) + (k * navirA * navirB), 0, -1.0, 1.0);
                K_I_aBc.reset();

                // Waux[IJk](A,Bc) += \sum(m) l_Im^Bc <Jk|Am>         (9)
                // Waux[IJk](A,Bc) += \sum(m) K[Jk](A,m) L[I](m,Bc)
                Waux->contract(false, false, navirA, navirA * navirB, naoccB, K_IjAk, l2AB,
                        (j * naoccB * navirA * naoccB) + (k * navirA * naoccB), (i * naoccB * navirA * navirB), 1.0, 1.0);

                // Waux[IJk](A,Bc) -= \sum(m) l_Jm^Bc <Ik|Am>         (10)
                // Waux[IJk](A,Bc) -= \sum(m) K[Ik](A,m) L[J](m,Bc)
                Waux->contract(false, false, navirA, navirA * navirB, naoccB, K_IjAk, l2AB,
                        (i * naoccB * navirA * naoccB) + (k * navirA * naoccB), (j * naoccB * navirA * navirB), -1.0, 1.0);

                // Wl_AAB(AB,c) += Waux(AB,c)
                Wl_AAB->axpy(Waux, 1.0);
                // Waux(A,Bc) sort 213 --> Waux(BA,c)
                // Wl_AAB(AB,c) -= Waux(BA,c)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int B = 0; B < navirA; ++B) {
                        Wl_AAB->axpy((size_t)navirB, A * navirA * navirB + B * navirB, 1, Waux, B * navirA * navirB + A * navirB, 1, -1.0);
                    }
                }
                Waux.reset();

                double Dijk = Dij + FockB->get(k + nfrzc, k + nfrzc);

                double Wijkabc, Wl_ijkabc, Vijkabc;
#pragma omp parallel for private(Wijkabc, Wl_ijkabc, Vijkabc) reduction(+ : sumAAB)
                for (long int a = 0; a < navirA; a++) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b < a; b++) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirB; c++) {
                            long int ac = ab_idxAB->get(a, c);
                            long int ca = ab_idxBA->get(c, a);
                            long int cb = ab_idxBA->get(c, b);
                            long int bc = ab_idxAB->get(b, c);

                            double Dijkabc = Dijkab - FockB->get(c + noccB, c + noccB);

                            // W[IJk](ABc) = W_AAB
                            //             + P_(AB) X4
                            Wijkabc = W_AAB->get(ab, c);
                            Wl_ijkabc = Wl_AAB->get(ab, c);

                            // V_AAB = t_Jk^Bc * f_IA + t_I^A * <Jk|Bc>
                            //       - t_Ik^Bc * f_JA - t_J^A * <Ik|Bc>
                            //       - t_Jk^Ac * f_IB - t_I^B * <Jk|Ac>
                            //       + t_Ik^Ac * f_JB + t_J^B * <Ik|Ac>
                            //       + t_IJ^AB * f_k^c + t_k^c * <IJ||AB>
                            //Vijkabc = t2AB->get(jk, bc) * FockA->get(nfrzc + i, noccA + a) + t1A->get(i, a) * K_IjAb->get(jk, bc)
                            //        - t2AB->get(ik, bc) * FockA->get(nfrzc + j, noccA + a) - t1A->get(j, a) * K_IjAb->get(ik, bc)
                            //        - t2AB->get(jk, ac) * FockA->get(nfrzc + i, noccA + b) - t1A->get(i, b) * K_IjAb->get(jk, ac)
                            //        + t2AB->get(ik, ac) * FockA->get(nfrzc + j, noccA + b) + t1A->get(j, b) * K_IjAb->get(ik, ac)
                            //        + t2AA->get(ij, ab) * FockB->get(nfrzc + k, noccB + c) + t1B->get(k, c) * G_IJAB->get(ij, ab);

                            // V_AAB = l_Jk^Bc * f_IA + l_I^A * <Jk|Bc>
                            //       - l_Ik^Bc * f_JA - l_J^A * <Ik|Bc>
                            //       - l_Jk^Ac * f_IB - l_I^B * <Jk|Ac>
                            //       + l_Ik^Ac * f_JB + l_J^B * <Ik|Ac>
                            //       + l_IJ^AB * f_k^c + l_k^c * <IJ||AB>
                            Vijkabc = l2AB->get(jk, bc) * FockA->get(nfrzc + i, noccA + a) + l1A->get(i, a) * K_IjAb->get(jk, bc)
                                    - l2AB->get(ik, bc) * FockA->get(nfrzc + j, noccA + a) - l1A->get(j, a) * K_IjAb->get(ik, bc)
                                    - l2AB->get(jk, ac) * FockA->get(nfrzc + i, noccA + b) - l1A->get(i, b) * K_IjAb->get(jk, ac)
                                    + l2AB->get(ik, ac) * FockA->get(nfrzc + j, noccA + b) + l1A->get(j, b) * K_IjAb->get(ik, ac)
                                    + l2AA->get(ij, ab) * FockB->get(nfrzc + k, noccB + c) + l1B->get(k, c) * G_IJAB->get(ij, ab);

                            // N9 : energy of AAB
                            //sumAAB += (Wijkabc + Vijkabc) * Wl_ijkabc / Dijkabc;
                            sumAAB += (Wl_ijkabc + Vijkabc) * Wijkabc / Dijkabc;

                        } // c
                    } // b
                } // a

                // progress counter
                ind += 1;
                double percent = static_cast<double>(ind) / static_cast<double>(Nijk) * 100.0;
                if (percent >= next_print) {
                    stop = std::time(nullptr);
                    next_print += step_print;
                    outfile->Printf("              %5.1lf  %8d s\n", percent,
                                    static_cast<int>(stop) - static_cast<int>(start));
                }
            } // k
        } // j
    } // i

    outfile->Printf("\tAAB (AT) energy                    : % 20.14f\n", sumAAB);
    E_at += sumAAB;

    // reset all AAB things
    G_IJAK.reset();
    G_IJAB.reset();
    W_AAB.reset();
    Wl_AAB.reset();
    t2AA.reset();
    l2AA.reset();
    K_I_aBc.reset();
    K_IjAb.reset();
    timer_off("CCSD(AT)-HM-AAB");

    timer_on("CCSD(AT)-HM-ABB");
    //==================================================
    //======================= ABB ======================
    //==================================================
    // N10 : beginning of ABB

    // malloc W_ABB[Ijk](Abc)
    SharedTensor2d W_ABB = std::make_shared<Tensor2d>("W[Ijk](A,cb)", navirA, navirB, navirB);

    // malloc Wl_ABB[Ijk](Abc)
    SharedTensor2d Wl_ABB = std::make_shared<Tensor2d>("Wl[Ijk](A,cb)", navirA, navirB, navirB);

    // read T2BB
    t2BB = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    t2BB->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // read L2BB
    l2BB = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    l2BB->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // progress counter
    start = std::time(nullptr);
    stop = std::time(nullptr);
    ind = 0;
    next_print = step_print;

    Nijk = naoccA * naoccB * (naoccB - 1) / 2;
    outfile->Printf("\n\tNumber of ijk combinations for ABB: %i \n", Nijk);

    // N11 : main loop of ABB
    double sumABB = 0.0;
    for (long int i = 0; i < naoccA; i++) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);

        for(long int j = 0; j < naoccB; j++) {
            long int ij = ij_idxAB->get(i, j);
            double Dij = Di + FockB->get(j + nfrzc, j + nfrzc);

            for (long int k = 0; k < j; k++) {
                long int kj = ij_idxBB->get(k, j);
                long int ik = ij_idxAB->get(i, k);
                long int jk = ij_idxBB->get(j, k);

                // read G[k](a,bc)
                G_i_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);
                G_i_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(k * navirB * navir2BB) * sizeof(double));

                // W_ABB[Ijk](A,cb) = \sum(e) t_Ij^Ae <ke||cb>          (1)
                // W_ABB[Ijk](A,cb) = T[Ij](A,e) G[k](e,cb)
                W_ABB->contract(false, false, navirA, navir2BB, navirB, t2AB, G_i_abc,
                        (i * naoccB * navirA * navirB) + (j * navirA * navirB), 0, 1.0, 0.0);

                // read G[j](a,bc)
                G_i_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(j * navirB * navir2BB) * sizeof(double));
                // W_ABB[Ijk](A,cb) = \sum(e) t_Ik^Ae <je||cb>          (2)
                // W_ABB[Ijk](A,cb) = T[Ik](A,e) G[j](e,cb)
                W_ABB->contract(false, false, navirA, navir2BB, navirB, t2AB, G_i_abc,
                        (i * naoccB * navirA * navirB) + (k * navirA * navirB), 0, -1.0, 1.0);
                G_i_abc.reset();

                // W_ABB[Ijk](A,cb) = \sum(m) t_jm^cb <Ik|Am>           (3)
                // W_ABB[Ijk](A,cb) = K[Ik](A,m) T[j](m,cb)
                W_ABB->contract(false, false, navirA, navir2BB, naoccB, K_IjAk, t2BB,
                        (i * naoccB * navirA * naoccB) + (k * navirA * naoccB), (j * naoccB * navir2BB), 1.0, 1.0);

                // W_ABB[Ijk](A,cb) = \sum(m) t_km^cb <Ij|Am>           (4)
                // W_ABB[Ijk](A,cb) = K[Ij](A,m) T[k](m,cb)
                W_ABB->contract(false, false, navirA, navir2BB, naoccB, K_IjAk, t2BB,
                        (i * naoccB * navirA * naoccB) + (j * navirA * naoccB), (k * naoccB * navir2BB), -1.0, 1.0);


                Waux = std::make_shared<Tensor2d>("X[Ijk](a,Bc)", navirB, navirA, navirB);

                // Waux[Ijk](b,Ac) = \sum(m) t_Im^Ac <kj||bm>       (5)
                // Waux[Ijk](b,Ac) = G[kj](b,m) T[I](m,Ac)
                Waux->contract(false, false, navirB, navirA * navirB, naoccB, G_ijak, t2AB,
                        (k * naoccB * navirB * naoccB) + (j * navirB * naoccB), (i * naoccB * navirA * navirB), 1.0, 0.0);

                // Waux(b,Ac) sort 231 --> Waux(A,cb)
                // Waux(b,Ac) sort 213 --> Waux(A,bc)
                // W_ABB(A,cb) += Waux(A,cb) - Waux(A,bc)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int b = 0; b < navirB; ++b) {
                        W_ABB->axpy((size_t)navirB, b * navirA * navirB + A * navirB, 1, Waux, A * navirB * navirB + b, navirB, 1.0);
                        W_ABB->axpy((size_t)navirB, b * navirA * navirB + A * navirB, 1, Waux, A * navirB * navirB + b * navirB, 1, -1.0);
                    }
                }

                // read K[I](a,Bc)
                K_I_aBc = std::make_shared<Tensor2d>("K[I] <a|Bc>", navirB, navirA, navirB);
                K_I_aBc->myread(psio_, PSIF_DFOCC_IABC_ABAB, (size_t)(i * navirB * navirA * navirB) * sizeof(double));

                // Waux[Ijk](c,Ab) = \sum(e) t_kj^ce <Ie|Ab>        (6)
                // Waux[Ijk](c,Ab) = T[kj](c,e) K[I](e,Ab)
                Waux->contract(false, false, navirB, navirA * navirB, navirB, t2BB, K_I_aBc,
                        (k * naoccB * navir2BB) + (j * navir2BB), 0, 1.0, 0.0);
                K_I_aBc.reset();

                // Waux(c,Ab) sort 213 --> Waux(A,cb)
                // Waux(c,Ab) sort 231 --> Waux(A,bc)
                // W_ABB(A,cb) += Waux(A,cb) - Waux(A,bc)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int c = 0; c < navirB; ++c) {
                        W_ABB->axpy((size_t)navirB, c * navirA * navirB + A * navirB, 1, Waux, A * navirB * navirB + c * navirB, 1, 1.0);
                        W_ABB->axpy((size_t)navirB, c * navirA * navirB + A * navirB, 1, Waux, A * navirB * navirB + c, navirB, -1.0);
                    }
                }
                Waux.reset();

                Waux = std::make_shared<Tensor2d>("X[Ijk](a,bC)", navirB, navirB, navirA);

                // read K[k](A,bC)
                K_i_AbC = std::make_shared<Tensor2d>("K[i] <A|bC>", navirA, navirB, navirA);
                K_i_AbC->myread(psio_, PSIF_DFOCC_IABC_BABA, (size_t)(k * navirB * navir2AA) * sizeof(double));

                // Waux[Ijk](b,cA) = \sum(E) t_Ij^Eb <kE|cA>        (7)
                // Waux[Ijk](b,cA) = T[Ij](E,b) K[k](E,cA)
                Waux->contract(true, false, navirB, navirB * navirA, navirA, t2AB, K_i_AbC,
                        (i * naoccB * navirA * navirB) + (j * navirA * navirB), 0, 1.0, 0.0);

                // read K[j](A,bC)
                K_i_AbC->myread(psio_, PSIF_DFOCC_IABC_BABA, (size_t)(j * navirB * navir2AA) * sizeof(double));

                // Waux[Ijk](b,cA) = \sum(E) t_Ik^Eb <jE|cA>        (8)
                // Waux[Ijk](b,cA) = T[Ik](E,b) K[j](E,cA)
                Waux->contract(true, false, navirB, navirB * navirA, navirA, t2AB, K_i_AbC,
                        (i * naoccB * navirA * navirB) + (k * navirA * navirB), 0, -1.0, 1.0);
                K_i_AbC.reset();

                // Waux(b,cA) sort 321 --> Waux(A,cb)
                // Waux(b,cA) sort 312 --> Waux(A,bc)
                // W_ABB(A,cb) += Waux(A,cb) - Waux(A,bc)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int b = 0; b < navirB; ++b) {
                        W_ABB->axpy((size_t)navirB, b * navirB * navirA + A, navirA, Waux, A * navirB * navirB + b, navirB, 1.0);
                        W_ABB->axpy((size_t)navirB, b * navirB * navirA + A, navirA, Waux, A * navirB * navirB + b * navirB, 1, -1.0);
                    }
                }

                // Waux[Ijk](c,bA) = \sum(M) t_kM^bA <Ij|Mc>        (9)
                // Waux[Ijk](c,bA) = K[Ij](M,c) t[k](M,bA)
                Waux->contract(true, false, navirB, navirB * navirA, naoccA, K_IjKa, t2BA,
                        (i * naoccB * naoccA * navirB) + (j * naoccA * navirB), (k * naoccA * navirB * navirA), 1.0, 0.0);

                // Waux[Ijk](c,bA) = \sum(M) t_jM^bA <Ik|Mc>        (10)
                // Waux[Ijk](c,bA) = K[Ik](M,c) t[j](M,bA)
                Waux->contract(true, false, navirB, navirB * navirA, naoccA, K_IjKa, t2BA,
                        (i * naoccB * naoccA * navirB) + (k * naoccA * navirB), (j * naoccA * navirB * navirA), -1.0, 1.0);

                // Waux(c,bA) sort 312 --> Waux(A,cb)
                // Waux(c,bA) sort 321 --> Waux(A,bc)
                // W_ABB(A,cb) += Waux(A,cb) - Waux(A,bc)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int c = 0; c < navirB; ++c) {
                        W_ABB->axpy((size_t)navirB, c * navirB * navirA + A, navirA, Waux, A * navirB * navirB + c * navirB, 1, 1.0);
                        W_ABB->axpy((size_t)navirB, c * navirB * navirA + A, navirA, Waux, A * navirB * navirB + c, navirB, -1.0);
                    }
                }
                Waux.reset();

                // ==========================
                // === Asymmetric Triples ===
                // ==========================

                // read G[k](a,bc)
                G_i_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);
                G_i_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(k * navirB * navir2BB) * sizeof(double));

                // Wl_ABB[Ijk](A,cb) = \sum(e) l_Ij^Ae <ke||cb>          (1)
                // Wl_ABB[Ijk](A,cb) = L[Ij](A,e) G[k](e,cb)
                Wl_ABB->contract(false, false, navirA, navir2BB, navirB, l2AB, G_i_abc,
                        (i * naoccB * navirA * navirB) + (j * navirA * navirB), 0, 1.0, 0.0);

                // read G[j](a,bc)
                G_i_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(j * navirB * navir2BB) * sizeof(double));
                // Wl_ABB[Ijk](A,cb) = \sum(e) l_Ik^Ae <je||cb>          (2)
                // Wl_ABB[Ijk](A,cb) = L[Ik](A,e) G[j](e,cb)
                Wl_ABB->contract(false, false, navirA, navir2BB, navirB, l2AB, G_i_abc,
                        (i * naoccB * navirA * navirB) + (k * navirA * navirB), 0, -1.0, 1.0);
                G_i_abc.reset();

                // Wl_ABB[Ijk](A,cb) = \sum(m) l_jm^cb <Ik|Am>           (3)
                // Wl_ABB[Ijk](A,cb) = K[Ik](A,m) L[j](m,cb)
                Wl_ABB->contract(false, false, navirA, navir2BB, naoccB, K_IjAk, l2BB,
                        (i * naoccB * navirA * naoccB) + (k * navirA * naoccB), (j * naoccB * navir2BB), 1.0, 1.0);

                // Wl_ABB[Ijk](A,cb) = \sum(m) l_km^cb <Ij|Am>           (4)
                // Wl_ABB[Ijk](A,cb) = K[Ij](A,m) L[k](m,cb)
                Wl_ABB->contract(false, false, navirA, navir2BB, naoccB, K_IjAk, l2BB,
                        (i * naoccB * navirA * naoccB) + (j * navirA * naoccB), (k * naoccB * navir2BB), -1.0, 1.0);


                Waux = std::make_shared<Tensor2d>("X[Ijk](a,Bc)", navirB, navirA, navirB);

                // Waux[Ijk](b,Ac) = \sum(m) l_Im^Ac <kj||bm>       (5)
                // Waux[Ijk](b,Ac) = G[kj](b,m) L[I](m,Ac)
                Waux->contract(false, false, navirB, navirA * navirB, naoccB, G_ijak, l2AB,
                        (k * naoccB * navirB * naoccB) + (j * navirB * naoccB), (i * naoccB * navirA * navirB), 1.0, 0.0);

                // Waux(b,Ac) sort 231 --> Waux(A,cb)
                // Waux(b,Ac) sort 213 --> Waux(A,bc)
                // Wl_ABB(A,cb) += Waux(A,cb) - Waux(A,bc)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int b = 0; b < navirB; ++b) {
                        Wl_ABB->axpy((size_t)navirB, b * navirA * navirB + A * navirB, 1, Waux, A * navirB * navirB + b, navirB, 1.0);
                        Wl_ABB->axpy((size_t)navirB, b * navirA * navirB + A * navirB, 1, Waux, A * navirB * navirB + b * navirB, 1, -1.0);
                    }
                }

                // read K[I](a,Bc)
                K_I_aBc = std::make_shared<Tensor2d>("K[I] <a|Bc>", navirB, navirA, navirB);
                K_I_aBc->myread(psio_, PSIF_DFOCC_IABC_ABAB, (size_t)(i * navirB * navirA * navirB) * sizeof(double));

                // Waux[Ijk](c,Ab) = \sum(e) l_kj^ce <Ie|Ab>        (6)
                // Waux[Ijk](c,Ab) = L[kj](c,e) K[I](e,Ab)
                Waux->contract(false, false, navirB, navirA * navirB, navirB, l2BB, K_I_aBc,
                        (k * naoccB * navir2BB) + (j * navir2BB), 0, 1.0, 0.0);
                K_I_aBc.reset();

                // Waux(c,Ab) sort 213 --> Waux(A,cb)
                // Waux(c,Ab) sort 231 --> Waux(A,bc)
                // Wl_ABB(A,cb) += Waux(A,cb) - Waux(A,bc)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int c = 0; c < navirB; ++c) {
                        Wl_ABB->axpy((size_t)navirB, c * navirA * navirB + A * navirB, 1, Waux, A * navirB * navirB + c * navirB, 1, 1.0);
                        Wl_ABB->axpy((size_t)navirB, c * navirA * navirB + A * navirB, 1, Waux, A * navirB * navirB + c, navirB, -1.0);
                    }
                }
                Waux.reset();

                Waux = std::make_shared<Tensor2d>("X[Ijk](a,bC)", navirB, navirB, navirA);

                // read K[k](A,bC)
                K_i_AbC = std::make_shared<Tensor2d>("K[i] <A|bC>", navirA, navirB, navirA);
                K_i_AbC->myread(psio_, PSIF_DFOCC_IABC_BABA, (size_t)(k * navirB * navir2AA) * sizeof(double));

                // Waux[Ijk](b,cA) = \sum(E) l_Ij^Eb <kE|cA>        (7)
                // Waux[Ijk](b,cA) = L[Ij](E,b) K[k](E,cA)
                Waux->contract(true, false, navirB, navirB * navirA, navirA, l2AB, K_i_AbC,
                        (i * naoccB * navirA * navirB) + (j * navirA * navirB), 0, 1.0, 0.0);

                // read K[j](A,bC)
                K_i_AbC->myread(psio_, PSIF_DFOCC_IABC_BABA, (size_t)(j * navirB * navir2AA) * sizeof(double));

                // Waux[Ijk](b,cA) = \sum(E) l_Ik^Eb <jE|cA>        (8)
                // Waux[Ijk](b,cA) = L[Ik](E,b) K[j](E,cA)
                Waux->contract(true, false, navirB, navirB * navirA, navirA, l2AB, K_i_AbC,
                        (i * naoccB * navirA * navirB) + (k * navirA * navirB), 0, -1.0, 1.0);
                K_i_AbC.reset();

                // Waux(b,cA) sort 321 --> Waux(A,cb)
                // Waux(b,cA) sort 312 --> Waux(A,bc)
                // Wl_ABB(A,cb) += Waux(A,cb) - Waux(A,bc)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int b = 0; b < navirB; ++b) {
                        Wl_ABB->axpy((size_t)navirB, b * navirB * navirA + A, navirA, Waux, A * navirB * navirB + b, navirB, 1.0);
                        Wl_ABB->axpy((size_t)navirB, b * navirB * navirA + A, navirA, Waux, A * navirB * navirB + b * navirB, 1, -1.0);
                    }
                }

                // Waux[Ijk](c,bA) = \sum(M) l_kM^bA <Ij|Mc>        (9)
                // Waux[Ijk](c,bA) = K[Ij](M,c) L[k](M,bA)
                Waux->contract(true, false, navirB, navirB * navirA, naoccA, K_IjKa, l2BA,
                        (i * naoccB * naoccA * navirB) + (j * naoccA * navirB), (k * naoccA * navirB * navirA), 1.0, 0.0);

                // Waux[Ijk](c,bA) = \sum(M) l_jM^bA <Ik|Mc>        (10)
                // Waux[Ijk](c,bA) = K[Ik](M,c) L[j](M,bA)
                Waux->contract(true, false, navirB, navirB * navirA, naoccA, K_IjKa, l2BA,
                        (i * naoccB * naoccA * navirB) + (k * naoccA * navirB), (j * naoccA * navirB * navirA), -1.0, 1.0);

                // Waux(c,bA) sort 312 --> Waux(A,cb)
                // Waux(c,bA) sort 321 --> Waux(A,bc)
                // Wl_ABB(A,cb) += Waux(A,cb) - Waux(A,bc)
//#pragma omp parallel for
                for (long int A = 0; A < navirA; ++A) {
                    for (long int c = 0; c < navirB; ++c) {
                        Wl_ABB->axpy((size_t)navirB, c * navirB * navirA + A, navirA, Waux, A * navirB * navirB + c * navirB, 1, 1.0);
                        Wl_ABB->axpy((size_t)navirB, c * navirB * navirA + A, navirA, Waux, A * navirB * navirB + c, navirB, -1.0);
                    }
                }
                Waux.reset();

                double Dijk = Dij + FockB->get(k + nfrzc, k + nfrzc);

                // read <ij||ab>
                G_ijab = std::make_shared<Tensor2d>("G <ij||ab>", naoccB, naoccB, navirB, navirB);
                G_ijab->read(psio_, PSIF_DFOCC_IJAB_BBBB);

                K_IjAb = std::make_shared<Tensor2d>("K <Ij|Ab>", naoccA, naoccB, navirA, navirB);
                K_IjAb->read(psio_, PSIF_DFOCC_IJAB_ABAB);

                double Wijkabc, Wl_ijkabc, Vijkabc;
#pragma omp parallel for private(Wijkabc, Wl_ijkabc, Vijkabc) reduction(+ : sumABB)
                for (long int a = 0; a < navirA; a++) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b < navirB; b++) {
                        double Dijkab = Dijka - FockB->get(b + noccB, b + noccB);
                        long int ab = ab_idxAB->get(a, b);
                        long int ba = ab_idxBA->get(b, a);
                        for (long int c = 0; c < b; c++) {
                            long int ac = ab_idxAB->get(a, c);
                            long int ca = ab_idxBA->get(c, a);
                            long int cb = ab_idxBB->get(c, b);
                            long int bc = ab_idxBB->get(b, c);

                            double Dijkabc = Dijkab - FockB->get(c + noccB, c + noccB);

                            // W_ABB[Ijk](Abc)
                            Wijkabc = W_ABB->get(a, cb);

                            // Wl_ABB[Ijk](Abc)
                            Wl_ijkabc = Wl_ABB->get(a, cb);

                            // V_ABB = t_Ij^Ab * f_kc + t_k^c * <Ij|Ab>
                            //       - t_Ik^Ab * f_jc - t_j^c * <Ik|Ab>
                            //       - t_Ij^Ac * f_kb - t_k^b * <Ij|Ac>
                            //       + t_Ik^Ac * f_jb + t_j^b * <Ik|Ac>
                            //       + t_kj^cb * f_I^A + t_I^A * <kj||cb>
                            //Vijkabc = t2AB->get(ij, ab) * FockB->get(nfrzc + k, noccB + c) + t1B->get(k, c) * K_IjAb->get(ij, ab)
                            //        - t2AB->get(ik, ab) * FockB->get(nfrzc + j, noccB + c) - t1B->get(j, c) * K_IjAb->get(ik, ab)
                            //        - t2AB->get(ij, ac) * FockB->get(nfrzc + k, noccB + b) - t1B->get(k, b) * K_IjAb->get(ij, ac)
                            //        + t2AB->get(ik, ac) * FockB->get(nfrzc + j, noccB + b) + t1B->get(j, b) * K_IjAb->get(ik, ac)
                            //        + t2BB->get(kj, cb) * FockA->get(nfrzc + i, noccA + a) + t1A->get(i, a) * G_ijab->get(kj, cb);

                            // V_ABB = l_Ij^Ab * f_kc + l_k^c * <Ij|Ab>
                            //       - l_Ik^Ab * f_jc - l_j^c * <Ik|Ab>
                            //       - l_Ij^Ac * f_kb - l_k^b * <Ij|Ac>
                            //       + l_Ik^Ac * f_jb + l_j^b * <Ik|Ac>
                            //       + l_kj^cb * f_I^A + l_I^A * <kj||cb>
                            Vijkabc = l2AB->get(ij, ab) * FockB->get(nfrzc + k, noccB + c) + l1B->get(k, c) * K_IjAb->get(ij, ab)
                                    - l2AB->get(ik, ab) * FockB->get(nfrzc + j, noccB + c) - l1B->get(j, c) * K_IjAb->get(ik, ab)
                                    - l2AB->get(ij, ac) * FockB->get(nfrzc + k, noccB + b) - l1B->get(k, b) * K_IjAb->get(ij, ac)
                                    + l2AB->get(ik, ac) * FockB->get(nfrzc + j, noccB + b) + l1B->get(j, b) * K_IjAb->get(ik, ac)
                                    + l2BB->get(kj, cb) * FockA->get(nfrzc + i, noccA + a) + l1A->get(i, a) * G_ijab->get(kj, cb);

                            // N12 : energy of ABB
                            //sumABB += (Wijkabc + Vijkabc) * Wl_ijkabc / Dijkabc;
                            sumABB += (Wl_ijkabc + Vijkabc) * Wijkabc / Dijkabc;

                        } // c
                    } // b
                } // a
                G_ijab.reset();
                K_IjAb.reset();

                // progress counter
                ind += 1;
                double percent = static_cast<double>(ind) / static_cast<double>(Nijk) * 100.0;
                if (percent >= next_print) {
                    stop = std::time(nullptr);
                    next_print += step_print;
                    outfile->Printf("              %5.1lf  %8d s\n", percent,
                                    static_cast<int>(stop) - static_cast<int>(start));
                }
            } // k
        } // j
    } // i


    outfile->Printf("\tABB (AT) energy                    : % 20.14f\n\n", sumABB);
    E_at += sumABB;

    // set energy
    Eccsd_at = Eccsd + E_at;

    // reset all ABB things
    W_ABB.reset();
    Wl_ABB.reset();
    G_ijak.reset();
    K_IjAk.reset();
    K_IjKa.reset();
    t2AB.reset();
    t2BA.reset();
    t2BB.reset();
    l2AB.reset();
    l2BA.reset();
    l2BB.reset();
    timer_off("CCSD(AT)-HM-ABB");

    // remove files
    remove_binary_file(PSIF_DFOCC_IABC_AAAA);
    remove_binary_file(PSIF_DFOCC_IABC_BBBB);
    remove_binary_file(PSIF_DFOCC_IABC_BABA);
    remove_binary_file(PSIF_DFOCC_IABC_ABAB);
} // uccsd_triples_hm()

}  // namespace dfoccwave
}  // namespace psi
