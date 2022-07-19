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

void DFOCC::uccsd_triples_grad_hm()
{
    pair_index();

    outfile->Printf("\tUsing high-memory disk algorithm...\n\n");

    timer_on("CCSD(T)-HM-AAA");
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

    // malloc G[I](ABC)
    G_I_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);
    SharedTensor2d G_J_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);

    // malloc tL2AA
    SharedTensor2d tL2AA = std::make_shared<Tensor2d>("(T)L2 <IJ||AB>", naoccA, naoccA, navirA, navirA);
    SharedTensor2d tL1A = std::make_shared<Tensor2d>("(T)L <I|A>", naoccA, navirA);

    // malloc M intr.
    SharedTensor2d M_OOVO = std::make_shared<Tensor2d>("M <IJ|AM>", naoccA, naoccA, navirA, naoccA);
    SharedTensor2d Mterm40s = std::make_shared<Tensor2d>("Mterm40s", navirA, navirA, navirA);

    // malloc correlation OPDM
    G1c_iiA = std::make_shared<Tensor1d>("(T) Correlation OPDM <I|I>", naoccA);
    G1c_aaA = std::make_shared<Tensor1d>("(T) Correlation OPDM <A|A>", navirA);

    // read T2AA
    SharedTensor2d t2AA = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    t2AA->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // progress counter
    std::time_t stop, start = std::time(nullptr);
    long int ind = 0;
    double step_print = 10.0;
    double next_print = step_print;

    long int Nijk = naoccA * naoccA * naoccA;
    outfile->Printf("\tNumber of ijk combinations for AAA: %i \n", Nijk);

    // N2 : main loop of AAA
    E_t = 0.0;
    double sumAAA = 0.0;
    for (long int i = 0; i < naoccA; i++) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);

        // read G[I](A,BC)
        G_I_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(i * navirA * navir2AA) * sizeof(double));
        Mterm40s->zero();
        for(long int j = 0; j < naoccA; j++) {
            long int ij = ij_idxAA->get(i, j);
            double Dij = Di + FockA->get(j + nfrzc, j + nfrzc);

            // read G[J](A,BC)
            G_J_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(j * navirA * navir2AA) * sizeof(double));
            for (long int k = 0; k < naoccA; k++) {
                long int kj = ij_idxAA->get(k, j);
                long int jk = ij_idxAA->get(j, k);
                long int ik = ij_idxAA->get(i, k);

                // malloc W[IJK](ABC)
                SharedTensor2d X_AAA = std::make_shared<Tensor2d>("X_AAA[IJK](A,CB)", navirA, navirA, navirA);

                // X_AAA[IJK](A,CB) -= \sum(E) t_KJ^AE <IE||CB>         (2)
                // X_AAA[IJK](A,CB) -= \sum(E) T[KJ](A,E) G[I](E,CB)
                X_AAA->contract(false, false, navirA, navir2AA, navirA, t2AA, G_I_ABC,
                        (k * naoccA * navir2AA) + (j * navir2AA), 0, -1.0, 0.0);

                // X_AAA[IJK](A,CB) -= \sum(E) t_IK^AE <JE||CB>         (3)
                // X_AAA[IJK](A,CB) -= \sum(E) T[IK](A,E) G[J](E,CB)
                X_AAA->contract(false, false, navirA, navir2AA, navirA, t2AA, G_J_ABC,
                        (i * naoccA * navir2AA) + (k * navir2AA), 0, -1.0, 1.0);

                // read G[K](A,BC)
                SharedTensor2d G_K_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);
                G_K_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(k * navirA * navir2AA) * sizeof(double));

                // X_AAA[IJK](A,CB) = \sum(E) t_IJ^AE <KE||CB>          (1)
                // X_AAA[IJK](A,CB) = \sum(E) T[IJ](A,E) G[K](E,CB)
                X_AAA->contract(false, false, navirA, navir2AA, navirA, t2AA, G_K_ABC,
                        (i * naoccA * navir2AA) + (j * navir2AA), 0, 1.0, 1.0);
                G_K_ABC.reset();

                SharedTensor2d Y_AAA = std::make_shared<Tensor2d>("Y_AAA[IJK](C,AB)", navirA, navirA, navirA);

                // Y_AAA[IJK](C,AB) = - \sum(M) t_IM^AB <KJ||CM>        (4)
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

                double Dijk = Dij + FockA->get(k + nfrzc, k + nfrzc);

                // malloc Tc[IJK](ABC)
                SharedTensor2d Tc_AAA = std::make_shared<Tensor2d>("Tc[IJK](A,BC)", navirA, navirA, navirA);
                // malloc Td[IJK](ABC)
                SharedTensor2d Td_AAA = std::make_shared<Tensor2d>("Td[IJK](A,BC)", navirA, navirA, navirA);

                double Wijkabc, Vijkabc;
//#pragma omp parallel for private(Wijkabc, Vijkabc) reduction(+ : sumAAA)
                for (long int a = 0; a < navirA; a++) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    double LIA = 0;
                    for (long int b = 0; b < navirA; b++) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; c++) {
                            long int ac = ab_idxAA->get(a, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            long int bc = ab_idxAA->get(b, c);

                            double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);

                            // W[IJK](ABC) = P(A/BC) X_AAA[IJK](ABC) + P(AB/C) Y_AAA[IJK](ABC)
                            Wijkabc = X_AAA->get(a, cb) - X_AAA->get(b, ca) - X_AAA->get(c, ab)
                                    + Y_AAA->get(c, ab) - Y_AAA->get(a, cb) - Y_AAA->get(b, ac);

                            // V_AAA = t_K^C <IJ||AB> + t_IJ^AB f_KC
                            //       - t_I^C <KJ||AB> - t_KJ^AB f_IC
                            //       - t_J^C <IK||AB> - t_IK^AB f_JC
                            //       - t_K^A <IJ||CB> + t_IJ^CB f_KA
                            //       + t_I^A <KJ||CB> - t_KJ^CB f_IA
                            //       + t_J^A <IK||CB> - t_IK^CB f_JA
                            //       - t_K^B <IJ||AC> + t_IJ^AC f_KB
                            //       + t_I^B <KJ||AC> - t_KJ^AC f_IB
                            //       + t_J^B <IK||AC> - t_IK^AC f_JB
                            Vijkabc = t1A->get(k, c) * G_IJAB->get(ij, ab) + t2AA->get(ij, ab) * FockA->get(nfrzc + k, noccA + c)
                                    - t1A->get(i, c) * G_IJAB->get(kj, ab) - t2AA->get(kj, ab) * FockA->get(nfrzc + i, noccA + c)
                                    - t1A->get(j, c) * G_IJAB->get(ik, ab) - t2AA->get(ik, ab) * FockA->get(nfrzc + j, noccA + c)
                                    - t1A->get(k, a) * G_IJAB->get(ij, cb) - t2AA->get(ij, cb) * FockA->get(nfrzc + k, noccA + a)
                                    + t1A->get(i, a) * G_IJAB->get(kj, cb) + t2AA->get(kj, cb) * FockA->get(nfrzc + i, noccA + a)
                                    + t1A->get(j, a) * G_IJAB->get(ik, cb) + t2AA->get(ik, cb) * FockA->get(nfrzc + j, noccA + a)
                                    - t1A->get(k, b) * G_IJAB->get(ij, ac) - t2AA->get(ij, ac) * FockA->get(nfrzc + k, noccA + b)
                                    + t1A->get(i, b) * G_IJAB->get(kj, ac) + t2AA->get(kj, ac) * FockA->get(nfrzc + i, noccA + b)
                                    + t1A->get(j, b) * G_IJAB->get(ik, ac) + t2AA->get(ik, ac) * FockA->get(nfrzc + j, noccA + b);

                            double tc_ijkabc = Wijkabc / Dijkabc;
                            Tc_AAA->set(a, bc, tc_ijkabc);

                            LIA += (0.25) * G_IJAB->get(jk, bc) * tc_ijkabc;

                            double td_ijkabc = Vijkabc / Dijkabc;
                            Td_AAA->set(a, bc, td_ijkabc);

                            // N3 : energy of AAA
                            sumAAA += (1.0/36.0) * (Wijkabc + Vijkabc) * Wijkabc / Dijkabc;

                        } // c
                    } // b
                    tL1A->add(i, a, LIA);
                } // a

                X_AAA.reset();
                Y_AAA.reset();

                // tijkabc = t(c)ijkabc + t(d)ijkabc
                Td_AAA->axpy(Tc_AAA, 1.0);

                // UCCSD(T) PDM Eq. (3) and (4)
                for (int a = 0; a < navirA; a++) {
                    double G1c_AAA_contrib = 0;
                    for (int b = 0; b < navirA; b++) {
                        for (int c = 0; c < navirA; c++) {
                            long int bc = ab_idxAA->get(b, c);
                            G1c_AAA_contrib += (1.0/12.0) * Td_AAA->get(a, bc) * Tc_AAA->get(a, bc);
                        }
                    }
                    G1c_iiA->subtract(i, G1c_AAA_contrib);
                    G1c_aaA->add(a, G1c_AAA_contrib);
                }

                // Mijkabc = tijkabc + t(c)ijkabc
                Td_AAA->axpy(Tc_AAA, 1.0);
                bool first_write = !(i || j || k);

                // UCCSD(T) M intr. Eq. (38)
                //------------------------------------------------------------------------------------------
                // Mterm38      = (1/2) *  \sum(BC) M[IJK](A,BC) * T(MK,BC)
                // Mterm38      = (1/2) * -\sum(BC) M[IJK](A,BC) * T(KM,BC)
                SharedTensor2d Mtemp = std::make_shared<Tensor2d>("Mtemp", navirA * navirA, navirA);
                Mtemp->trans(Td_AAA);
                // Mterm38(M,A) = (1/2) * -\sum(BC) T[K](M,BC)   * M[IJK](BC,A)
                SharedTensor2d Mterm38 = std::make_shared<Tensor2d>("Mterm38", naoccA, navirA);
                Mterm38->contract(false, false, naoccA, navirA, navir2AA, t2AA, Mtemp,
                        (k * naoccA * navir2AA), 0, -0.5, 0.0);
                Mtemp.reset();

                SharedTensor2d Mterm38tr = std::make_shared<Tensor2d>("Mterm38tr", navirA, naoccA);
                Mterm38tr->trans(Mterm38);
                Mterm38.reset();

                // M_OOVO(IJ,AM) += \sum(K) Mterm38tr(A,M)
                M_OOVO->axpy((size_t)(naoccA * navirA), 0, 1, Mterm38tr, (i * naoccA * navirA * naoccA) + (j * navirA * naoccA), 1, 1.0);
                Mterm38tr.reset();

                // UCCSD(T) M intr. Eq. (39)
                //------------------------------------------------------------------------------------------
                // Mterm39(AB) = \sum(C) Tc[IJK](AB,C) * T[K](C)
                SharedTensor1d Mterm39 = std::make_shared<Tensor1d>("Mterm39", navir2AA);
                Mterm39->gemv(false, navir2AA, navirA, Tc_AAA, t1A,
                        0, (k * navirA), 1.0, 0.0);

                // M_OOVV(IJ,AB) += \sum(K) Mterm39(AB)
                SharedTensor2d M_OOVV = std::make_shared<Tensor2d>("M <IJ||AB>", naoccA, naoccA, navirA, navirA);
                if (!first_write) M_OOVV->read_anti_symm(psio_, PSIF_DFOCC_DENS);
                M_OOVV->axpy((size_t)(navirA * navirA), 0, 1, Mterm39, (i * naoccA * navirA * navirA) + (j * navirA * navirA), 1, 1.0);
                M_OOVV->write_anti_symm(psio_, PSIF_DFOCC_DENS);
                M_OOVV.reset();
                Mterm39.reset();

                Tc_AAA.reset();

                // UCCSD(T) Lambda Eq. (20) 1st term
                //------------------------------------------------------------------------------------------
                // read G[K](A,BC)
                G_K_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);
                G_K_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(k * navirA * navir2AA) * sizeof(double));

                // lterm20(A,B) = (1/2) * -\sum(EF) M[IJK](A,EF) G[K](B,EF)
                // lterm20(A,B) = (1/2) * -\sum(EF) M[IJK](A,EF) G[K](EF,B)
                SharedTensor2d lterm20 = std::make_shared<Tensor2d>("lterm20", navirA, navirA);
                SharedTensor2d lterm20tr = std::make_shared<Tensor2d>("lterm20tr", navirA, navirA);
                lterm20->contract(false, true, navirA, navirA, navir2AA, Td_AAA, G_K_ABC, -0.5, 0.0);
                G_K_ABC.reset();

                // tL2AA[IJ](AB) += P_(AB) \sum(K) lterm20(AB)
                tL2AA->axpy((size_t)navir2AA, 0, 1, lterm20, (i * naoccA * navir2AA) + (j * navir2AA), 1, 1.0);
                lterm20tr->trans(lterm20);
                lterm20.reset();
                tL2AA->axpy((size_t)navir2AA, 0, 1, lterm20tr, (i * naoccA * navir2AA) + (j * navir2AA), 1, -1.0);
                lterm20tr.reset();

                // UCCSD(T) Lambda Eq. (21)
                //------------------------------------------------------------------------------------------
                SharedTensor2d lterm21 = std::make_shared<Tensor2d>("lterm21", navir2AA, naoccA);
                // lterm21(AB,L) = (1/2) * \sum(C) M[IJK](AB,C) <JK||CL>
                // lterm21(AB,L) = (1/2) * \sum(C) M[IJK](AB,C) G[JK](C,L)
                lterm21->contract(false, false, navir2AA, naoccA, navirA, Td_AAA, G_IJAK,
                        0, (j * naoccA * navirA * naoccA) + (k * navirA * naoccA), 0.5, 0.0);

                SharedTensor2d lterm21_t = std::make_shared<Tensor2d>("lterm21_t", naoccA, navir2AA);
                lterm21_t->trans(lterm21);
                lterm21.reset();

                // tL2AA(IL,AB) += P_(IL) lterm21_t[I](L,AB)
                tL2AA->axpy((size_t)naoccA * navir2AA, 0, 1, lterm21_t, i * naoccA * navir2AA, 1, 1.0);

                for (int l = 0; l < naoccA; l++) {
                    long int li = ij_idxAA->get(l, i);
                    for (int a = 0; a < navirA; a++) {
                        for (int b = 0; b < navirA; b++) {
                            long int ab = ab_idxAA->get(a, b);
                            double val = lterm21_t->get(l, ab);
                            tL2AA->subtract(li, ab, val);
                        }
                    }
                }
                lterm21_t.reset();

                // UCCSD(T) M intr. Eq. (40)
                //------------------------------------------------------------------------------------------
                // Mterm40       = (1/2) * \sum(C) M[IJK](A,BC) * T(JK,DC)
                // Mterm40       = (1/2) * \sum(C) M[IJK](C,AB) * T(JK,DC)
                Mtemp = std::make_shared<Tensor2d>("Mtemp", navirA, navirA, navirA);
                Mtemp->sort3a(312, navirA, navirA, navirA, Td_AAA, 1.0, 0.0);
                // Mterm40(D,AB) = (1/2) * \sum(C) T[JK](D,C)   * M[IJK](C,AB)
                SharedTensor2d Mterm40 = std::make_shared<Tensor2d>("Mterm40", navirA, navirA, navirA);
                Mterm40->contract(false, false, navirA, navir2AA, navirA, t2AA, Mtemp,
                        (j * naoccA * navir2AA) + (k * navir2AA), 0, 0.5, 0.0);
                Mtemp.reset();

                Mterm40s->sort3a(231, navirA, navirA, navirA, Mterm40, 1.0, 1.0);
                Mterm40.reset();

                Td_AAA.reset();

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
        Mterm40s->mywrite(psio_, PSIF_DFOCC_MIABC_AAAA, (bool)i);
        Mterm40s->zero();
    } // i

    Mterm40s.reset();

    outfile->Printf("\tAAA (T) energy                     : % 20.14f\n", sumAAA);
    E_t += sumAAA;

    // reset all AAA things
    t2AA.reset();
    G_I_ABC.reset();
    G_J_ABC.reset();
    G_IJAB.reset();
    timer_off("CCSD(T)-HM-AAA");

    timer_on("CCSD(T)-HM-BBB");
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

    // malloc tL2BB
    SharedTensor2d tL1B = std::make_shared<Tensor2d>("(T)L <i|a>", naoccB, navirB);
    SharedTensor2d tL2BB = std::make_shared<Tensor2d>("(T)L2 <ij||ab>", naoccB, naoccB, navirB, navirB);

    // malloc M intr.
    SharedTensor2d M_oovo = std::make_shared<Tensor2d>("M <ij|am>", naoccB, naoccB, navirB, naoccB);
    SharedTensor2d Mterm43s = std::make_shared<Tensor2d>("Mterm43s", navirB, navirB, navirB);

    // malloc correlation OPDM
    G1c_iiB = std::make_shared<Tensor1d>("(T) Correlation OPDM <i|i>", naoccB);
    G1c_aaB = std::make_shared<Tensor1d>("(T) Correlation OPDM <a|a>", navirB);

    // read T2BB
    SharedTensor2d t2BB = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    t2BB->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // malloc G[i](a,bc)
    G_i_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);
    SharedTensor2d G_j_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);

    // progress counter
    start = std::time(nullptr);
    stop = std::time(nullptr);
    ind = 0;
    next_print = step_print;

    Nijk = naoccB * naoccB * naoccB;
    outfile->Printf("\n\tNumber of ijk combinations for BBB: %i \n", Nijk);

    // N5 : main loop of BBB
    double sumBBB = 0.0;
    for (long int i = 0; i < naoccB; i++) {
        double Di = FockB->get(i + nfrzc, i + nfrzc);

        Mterm43s->zero();
        // read G[i](a,bc)
        G_i_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(i * navirB * navir2BB) * sizeof(double));
        for(long int j = 0; j < naoccB; j++) {
            long int ij = ij_idxBB->get(i, j);
            double Dij = Di + FockB->get(j + nfrzc, j + nfrzc);

            // read G[j](a,bc)
            G_j_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(j * navirB * navir2BB) * sizeof(double));
            for (long int k = 0; k < naoccB; k++) {
                long int kj = ij_idxBB->get(k, j);
                long int jk = ij_idxBB->get(j, k);
                long int ik = ij_idxBB->get(i, k);

                // read G[k](a,bc)
                SharedTensor2d G_k_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);
                G_k_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(k * navirB * navir2BB) * sizeof(double));

                // malloc W[ijk](abc)
                SharedTensor2d X_BBB = std::make_shared<Tensor2d>("X[ijk](a,cb)", navirB, navirB, navirB);

                // X_BBB[ijk](a,cb) = \sum(e) t_ij^ae <ke||cb>          (1)
                // X_BBB[ijk](a,cb) = \sum(e) T[ij](a,e) G[k](e,cb)
                X_BBB->contract(false, false, navirB, navir2BB, navirB, t2BB, G_k_abc,
                        (i * naoccB * navir2BB) + (j * navir2BB), 0, 1.0, 0.0);
                G_k_abc.reset();

                // X_BBB[ijk](a,cb) -= \sum(e) t_kj^ae <ie||cb>         (2)
                // X_BBB[ijk](a,cb) -= \sum(e) T[kj](a,e) G[i](e,cb)
                X_BBB->contract(false, false, navirB, navir2BB, navirB, t2BB, G_i_abc,
                        (k * naoccB * navir2BB) + (j * navir2BB), 0, -1.0, 1.0);

                // X_BBB[ijk](a,cb) -= \sum(e) t_ik^ae <je||cb>         (3)
                // X_BBB[ijk](a,cb) -= \sum(e) T[ik](a,e) G[j](e,cb)
                X_BBB->contract(false, false, navirB, navir2BB, navirB, t2BB, G_j_abc,
                        (i * naoccB * navir2BB) + (k * navir2BB), 0, -1.0, 1.0);

                SharedTensor2d Y_BBB = std::make_shared<Tensor2d>("Y[ijk](c,ab)", navirB, navirB, navirB);

                // Y_BBB[ijk](c,ab) = - \sum(m) t_im^ab <kj||cm         (4)
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

                double Dijk = Dij + FockB->get(k + nfrzc, k + nfrzc);

                // malloc Tc[ijk](abc)
                SharedTensor2d Tc_BBB = std::make_shared<Tensor2d>("Tc[ijk](a,bc)", navirB, navirB, navirB);
                // malloc Td[ijk](abc)
                SharedTensor2d Td_BBB = std::make_shared<Tensor2d>("Td[ijk](a,bc)", navirB, navirB, navirB);

                double Wijkabc, Vijkabc;
//#pragma omp parallel for private(Wijkabc, Vijkabc) reduction(+ : sumBBB)
                for (long int a = 0; a < navirB; a++) {
                    double Dijka = Dijk - FockB->get(a + noccB, a + noccB);
                    double Lia = 0;
                    for (long int b = 0; b < navirB; b++) {
                        double Dijkab = Dijka - FockB->get(b + noccB, b + noccB);
                        long int ab = ab_idxBB->get(a, b);
                        long int ba = ab_idxBB->get(b, a);
                        for (long int c = 0; c < navirB; c++) {
                            long int ac = ab_idxBB->get(a, c);
                            long int ca = ab_idxBB->get(c, a);
                            long int cb = ab_idxBB->get(c, b);
                            long int bc = ab_idxBB->get(b, c);

                            // W[ijk](abc) = P(a/bc) X_BBB[ijk](abc) + P(ab/c) Y_BBB[ijk](abc)
                            Wijkabc = X_BBB->get(a, cb) - X_BBB->get(b, ca) - X_BBB->get(c, ab)
                                    + Y_BBB->get(c, ab) - Y_BBB->get(a, cb) - Y_BBB->get(b, ac);

                            // V_BBB = t_k^c <ij||ab> + t_ij^ab f_kc
                            //       - t_i^c <kj||ab> - t_kj^ab f_ic
                            //       - t_j^c <ik||ab> - t_ik^ab f_jc
                            //       - t_k^a <ij||cb> + t_ij^cb f_ka
                            //       + t_i^a <kj||cb> - t_kj^cb f_ia
                            //       + t_j^a <ik||cb> - t_ik^cb f_ja
                            //       - t_k^b <ij||ac> + t_ij^ac f_kb
                            //       + t_i^b <kj||ac> - t_kj^ac f_ib
                            //       + t_j^b <ik||ac> - t_ik^ac f_jb
                            Vijkabc = t1B->get(k, c) * G_ijab->get(ij, ab) + t2BB->get(ij, ab) * FockB->get(nfrzc + k, noccB + c)
                                    - t1B->get(i, c) * G_ijab->get(kj, ab) - t2BB->get(kj, ab) * FockB->get(nfrzc + i, noccB + c)
                                    - t1B->get(j, c) * G_ijab->get(ik, ab) - t2BB->get(ik, ab) * FockB->get(nfrzc + j, noccB + c)
                                    - t1B->get(k, a) * G_ijab->get(ij, cb) - t2BB->get(ij, cb) * FockB->get(nfrzc + k, noccB + a)
                                    + t1B->get(i, a) * G_ijab->get(kj, cb) + t2BB->get(kj, cb) * FockB->get(nfrzc + i, noccB + a)
                                    + t1B->get(j, a) * G_ijab->get(ik, cb) + t2BB->get(ik, cb) * FockB->get(nfrzc + j, noccB + a)
                                    - t1B->get(k, b) * G_ijab->get(ij, ac) - t2BB->get(ij, ac) * FockB->get(nfrzc + k, noccB + b)
                                    + t1B->get(i, b) * G_ijab->get(kj, ac) + t2BB->get(kj, ac) * FockB->get(nfrzc + i, noccB + b)
                                    + t1B->get(j, b) * G_ijab->get(ik, ac) + t2BB->get(ik, ac) * FockB->get(nfrzc + j, noccB + b);

                            double Dijkabc = Dijkab - FockB->get(c + noccB, c + noccB);

                            double tc_ijkabc = Wijkabc / Dijkabc;
                            Tc_BBB->set(a, bc, tc_ijkabc);

                            Lia += (0.25) * G_ijab->get(jk, bc) * tc_ijkabc;

                            double td_ijkabc = Vijkabc / Dijkabc;
                            Td_BBB->set(a, bc, td_ijkabc);

                            // N6 : energy of BBB
                            sumBBB += (1.0/36.0) * (Wijkabc + Vijkabc) * Wijkabc / Dijkabc;

                        } // c
                    } // b
                    tL1B->add(i, a, Lia);
                } // a

                X_BBB.reset();
                Y_BBB.reset();

                // tijkabc = t(c)ijkabc + t(d)ijkabc
                Td_BBB->axpy(Tc_BBB, 1.0);

                // UCCSD(T) PDM Eq. (5) and (8)
                for (int a = 0; a < navirB; a++) {
                    double G1c_BBB_contrib = 0;
                    for (int b = 0; b < navirB; b++) {
                        for (int c = 0; c < navirB; c++) {
                            long int bc = ab_idxBB->get(b, c);
                            G1c_BBB_contrib += (1.0/12.0) * Td_BBB->get(a, bc) * Tc_BBB->get(a, bc);
                        }
                    }
                    G1c_iiB->subtract(i, G1c_BBB_contrib);
                    G1c_aaB->add(a, G1c_BBB_contrib);
                }

                // Mijkabc = tijkabc + t(c)ijkabc
                Td_BBB->axpy(Tc_BBB, 1.0);
                bool first_write = !(i || j || k);

                // UCCSD(T) M intr. Eq. (41)
                //------------------------------------------------------------------------------------------
                // Mterm41      = (1/2) *  \sum(bc) M[ijk](a,bc) * T(mk,bc)
                // Mterm41      = (1/2) * -\sum(bc) M[ijk](a,bc) * T(km,bc)
                SharedTensor2d Mtemp = std::make_shared<Tensor2d>("Mtemp", navirB * navirB, navirB);
                Mtemp->trans(Td_BBB);
                // Mterm41(m,a) = (1/2) * -\sum(bc) T[k](m,bc)   * M[ijk](bc,a)
                SharedTensor2d Mterm41 = std::make_shared<Tensor2d>("Mterm41", naoccB, navirB);
                Mterm41->contract(false, false, naoccB, navirB, navir2BB, t2BB, Mtemp,
                        (k * naoccB * navir2BB), 0, -0.5, 0.0);
                Mtemp.reset();

                SharedTensor2d Mterm41tr = std::make_shared<Tensor2d>("Mterm41tr", navirB, naoccB);
                Mterm41tr->trans(Mterm41);
                Mterm41.reset();

                // M_oovo(ij,am) += \sum(k) Mterm41tr(a,m)
                M_oovo->axpy((size_t)(naoccB * navirB), 0, 1, Mterm41tr, (i * naoccB * navirB * naoccB) + (j * navirB * naoccB), 1, 1.0);
                Mterm41tr.reset();

                // UCCSD(T) M intr. Eq. (42)
                //------------------------------------------------------------------------------------------
                // Mterm42(ab) = \sum(c) Tc[ijk](ab,c) * T[k](c)
                // Checked (ijab)   <-->   (crawford T3_grad_UHF_BBB.cc Gijab) | 7
                SharedTensor1d Mterm42 = std::make_shared<Tensor1d>("Mterm42", navir2BB);
                Mterm42->gemv(false, navir2BB, navirB, Tc_BBB, t1B,
                        0, (k * navirB), 1.0, 0.0);

                // M_oovv(ij,ab) += \sum(k) Mterm42(ab)
                SharedTensor2d M_oovv = std::make_shared<Tensor2d>("M <ij||ab>", naoccB, naoccB, navirB, navirB);
                if (!first_write) M_oovv->read_anti_symm(psio_, PSIF_DFOCC_DENS);
                M_oovv->axpy((size_t)(navirB * navirB), 0, 1, Mterm42, (i * naoccB * navirB * navirB) + (j * navirB * navirB), 1, 1.0);
                M_oovv->write_anti_symm(psio_, PSIF_DFOCC_DENS);
                M_oovv.reset();
                Mterm42.reset();

                Tc_BBB.reset();

                // UCCSD(T) Lambda Eq. (22) 1st term
                //------------------------------------------------------------------------------------------
                G_k_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);
                G_k_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(k * navirB * navir2BB) * sizeof(double));

                // lterm22(a,b) = (1/2) * -\sum(ef) M[ijk](a,ef) <kb||ef>
                // lterm22(a,b) = (1/2) * -\sum(ef) M[ijk](a,ef) G[k](ef,b)
                SharedTensor2d lterm22 = std::make_shared<Tensor2d>("lterm22", navirB, navirB);
                SharedTensor2d lterm22_tr = std::make_shared<Tensor2d>("lterm22tr", navirB, navirB);
                lterm22->contract(false, true, navirB, navirB, navir2BB, Td_BBB, G_k_abc, -0.5, 0.0);
                G_k_abc.reset();

                // tL2BB[ij](ab) += P_(ab) \sum(k) lterm22(ab)
                tL2BB->axpy((size_t)navirB * navirB, 0, 1, lterm22, i * naoccB * navir2BB + j * navir2BB, 1, 1.0);
                lterm22_tr->trans(lterm22);
                lterm22.reset();
                tL2BB->axpy((size_t)navirB * navirB, 0, 1, lterm22_tr, i * naoccB * navir2BB + j * navir2BB, 1, -1.0);
                lterm22_tr.reset();

                // UCCSD(T) Lambda Eq. (23)
                //------------------------------------------------------------------------------------------
                SharedTensor2d lterm23 = std::make_shared<Tensor2d>("lterm23", navir2BB, naoccB);
                // lterm23(ab,l) = (1/2) * \sum(c) M[ijk](ab,c) <jk||cl>
                // lterm23(ab,l) = (1/2) * \sum(c) M[ijk](ab,c) G[jk](c,l)
                lterm23->contract(false, false, navir2BB, naoccB, navirB, Td_BBB, G_ijak,
                        0, (j * naoccB * navirB * naoccB) + (k * navirB * naoccB), 0.5, 0.0);

                SharedTensor2d lterm23_t = std::make_shared<Tensor2d>("lterm23_t", naoccB, navir2BB);
                lterm23_t->trans(lterm23);
                lterm23.reset();

                // tL2BB(il,ab) += \sum(jk) lterm23_t[ijk](l,ab)
                tL2BB->axpy((size_t)naoccB * navir2BB, 0, 1, lterm23_t, i * naoccB * navir2BB, 1, 1.0);

                for (int l = 0; l < naoccB; l++) {
                    long int li = ij_idxBB->get(l, i);
                    for (int a = 0; a < navirB; a++) {
                        for (int b = 0; b < navirB; b++) {
                            long int ab = ab_idxBB->get(a, b);
                            double val = lterm23_t->get(l, ab);
                            tL2BB->subtract(li, ab, val);
                        }
                    }
                }
                lterm23_t.reset();

                // UCCSD(T) M intr. Eq. (43)
                //------------------------------------------------------------------------------------------
                // Mterm43       = (1/2) * \sum(c) M[ijk](a,bc) * T(jk,dc)
                // Mterm43       = (1/2) * \sum(c) M[ijk](c,ab) * T(jk,dc)
                Mtemp = std::make_shared<Tensor2d>("Mtemp", navirB, navirB, navirB);
                Mtemp->sort3a(312, navirB, navirB, navirB, Td_BBB, 1.0, 0.0);
                // Mterm43(d,ab) = (1/2) * \sum(c) T[jk](d,c)   * M[ijk](c,ab)
                SharedTensor2d Mterm43 = std::make_shared<Tensor2d>("Mterm43", navirB, navirB, navirB);
                Mterm43->contract(false, false, navirB, navir2BB, navirB, t2BB, Mtemp,
                        (j * naoccB * navir2BB) + (k * navir2BB), 0, 0.5, 0.0);
                Mtemp.reset();

                Mterm43s->sort3a(231, navirB, navirB, navirB, Mterm43, 1.0, 1.0);
                Mterm43.reset();

                Td_BBB.reset();

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
        Mterm43s->mywrite(psio_, PSIF_DFOCC_MIABC_BBBB, (bool)i);
        Mterm43s->zero();
    } // i

    Mterm43s.reset();

    outfile->Printf("\tBBB (T) energy                     : % 20.14f\n", sumBBB);
    E_t += sumBBB;

    // reset all BBB things
    t2BB.reset();
    G_ijab.reset();
    G_i_abc.reset();
    timer_off("CCSD(T)-HM-BBB");

    //==================================================
    //======================= AAB ======================
    //==================================================
    // N7 : beginning of AAB
    timer_on("CCSD(T)-HM-AAB");

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
    SharedTensor2d K_IjaK = std::make_shared<Tensor2d>("K <Ij|aK>", naoccA, naoccB, navirB, naoccA);//UB
    K_IjaK->sort(1243, K_IjKa, 1.0, 0.0);//UB

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

    // malloc M intr.
    SharedTensor2d M_oOvO = std::make_shared<Tensor2d>("M <iJ|aM>", naoccB, naoccA, navirB, naoccA);
    SharedTensor2d M_OoVo = std::make_shared<Tensor2d>("M <Ij|Am>", naoccA, naoccB, navirA, naoccB);
    SharedTensor2d Mterm49s = std::make_shared<Tensor2d>("Mterm49s(A,cd)", navirA, navirB, navirB);
    SharedTensor2d Mterm45diskalloc = std::make_shared<Tensor2d>("Mterm45s(c,BD)", navirB, navirA, navirA);
    for (int i = 0; i < naoccB; i++) {
        Mterm45diskalloc->mywrite(psio_, PSIF_DFOCC_MIABC_BBAA, (bool)i);
    }

    // malloc tL2AB
    SharedTensor2d tL2AB = std::make_shared<Tensor2d>("(T)L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);

    G_IJAB = std::make_shared<Tensor2d>("G <IJ||AB>", naoccA, naoccA, navirA, navirA);
    G_IJAB->read(psio_, PSIF_DFOCC_IJAB_AAAA);

    // progress counter
    start = std::time(nullptr);
    stop = std::time(nullptr);
    ind = 0;
    next_print = step_print;

    Nijk = naoccA * naoccA * naoccB;
    outfile->Printf("\n\tNumber of ijk combinations for AAB: %i \n", Nijk);

    // N8 : main loop of AAB
    double sumAAB = 0.0;
    for (long int i = 0; i < naoccA; i++) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);

        for(long int j = 0; j < naoccA; j++) {
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

                // malloc Tc[ijk](abc)
                SharedTensor2d Tc_AAB = std::make_shared<Tensor2d>("Tc[IJk](A,Bc)", navirA, navirA, navirB);
                // malloc Td[ijk](abc)
                SharedTensor2d Td_AAB = std::make_shared<Tensor2d>("Td[IJk](A,Bc)", navirA, navirA, navirB);

                double Dijk = Dij + FockB->get(k + nfrzc, k + nfrzc);

                double Wijkabc, Vijkabc;
//#pragma omp parallel for private(Wijkabc, Vijkabc) reduction(+ : sumAAB)
                for (long int a = 0; a < navirA; a++) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    double LIA = 0;
                    for (long int b = 0; b < navirA; b++) {
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

                            // V_AAB = t_Jk^Bc * f_IA + t_I^A * <Jk|Bc>
                            //       - t_Ik^Bc * f_JA - t_J^A * <Ik|Bc>
                            //       - t_Jk^Ac * f_IB - t_I^B * <Jk|Ac>
                            //       + t_Ik^Ac * f_JB + t_J^B * <Ik|Ac>
                            //       + t_IJ^AB * f_k^c + t_k^c * <IJ||AB>
                            Vijkabc = t2AB->get(jk, bc) * FockA->get(nfrzc + i, noccA + a) + t1A->get(i, a) * K_IjAb->get(jk, bc)
                                    - t2AB->get(ik, bc) * FockA->get(nfrzc + j, noccA + a) - t1A->get(j, a) * K_IjAb->get(ik, bc)
                                    - t2AB->get(jk, ac) * FockA->get(nfrzc + i, noccA + b) - t1A->get(i, b) * K_IjAb->get(jk, ac)
                                    + t2AB->get(ik, ac) * FockA->get(nfrzc + j, noccA + b) + t1A->get(j, b) * K_IjAb->get(ik, ac)
                                    + t2AA->get(ij, ab) * FockB->get(nfrzc + k, noccB + c) + t1B->get(k, c) * G_IJAB->get(ij, ab);

                            double tc_ijkabc = Wijkabc / Dijkabc;
                            Tc_AAB->set(a, bc, tc_ijkabc);

                            LIA += K_IjAb->get(jk, bc) * tc_ijkabc;

                            double Lia = (0.25) * G_IJAB->get(ij, ab) * tc_ijkabc;
                            tL1B->add(k, c, Lia);

                            double td_ijkabc = Vijkabc / Dijkabc;
                            Td_AAB->set(a, bc, td_ijkabc);

                            // N9 : energy of AAB
                            sumAAB += (0.25) * (Wijkabc + Vijkabc) * Wijkabc / Dijkabc;

                        } // c

                    } // b
                    tL1A->add(i, a, LIA);
                } // a
                W_AAB.reset();

                // tijkabc = t(c)ijkabc + t(d)ijkabc
                Td_AAB->axpy(Tc_AAB, 1.0);

                // UCCSD(T) PDM Eq. (3) and (4)
                for (int a = 0; a < navirA; a++) {
                    double G1c_AAB_contrib = 0;
                    for (int b = 0; b < navirA; b++) {
                        for (int c = 0; c < navirB; c++) {
                            long int bc = ab_idxAB->get(b, c);
                            G1c_AAB_contrib += (1.0/2.0) * Td_AAB->get(a, bc) * Tc_AAB->get(a, bc);
                            double G1c_AAB_contrib2 = (1.0/4.0) * Td_AAB->get(a, bc) * Tc_AAB->get(a, bc);
                            G1c_iiB->subtract(k, G1c_AAB_contrib2);
                            G1c_aaB->add(c, G1c_AAB_contrib2);
                        }
                    }
                    G1c_iiA->subtract(i, G1c_AAB_contrib);
                    G1c_aaA->add(a, G1c_AAB_contrib);
                }

                // Mijkabc = tijkabc + t(c)ijkabc
                Td_AAB->axpy(Tc_AAB, 1.0);
                bool first_write = !(i || j || k);

                // UCCSD(T) M intr. Eq. (44)
                //------------------------------------------------------------------------------------------
                // Mterm44      = (1/2) *  \sum(BA) M[IJk](A,Bc) * T(MI,BA)
                // Mterm44      = (1/2) *  \sum(BA) M[IJk](A,Bc) * T(IM,AB)
                // Mterm44(M,c) = (1/2) *  \sum(BA) T[I]M,AB)    * M[IJk](AB,c) 
                SharedTensor2d Mterm44 = std::make_shared<Tensor2d>("Mterm44", naoccA, navirB);
                Mterm44->contract(false, false, naoccA, navirB, navir2AA, t2AA, Td_AAB,
                        (i * naoccA * navir2AA), 0, 0.5, 0.0);

                SharedTensor2d Mterm44tr = std::make_shared<Tensor2d>("Mterm44tr(c,M)", navirB, naoccA);
                Mterm44tr->trans(Mterm44);
                Mterm44.reset();

                // M_oOvO[kJ](cM) += \sum(I) Mterm44tr(cM)
                M_oOvO->axpy((size_t)(navirB * naoccA), 0, 1, Mterm44tr, (k * naoccA * navirB * naoccA) + (j * navirB * naoccA), 1, 1.0);
                Mterm44tr.reset();

                // UCCSD(T) M intr. Eq. (46)
                //------------------------------------------------------------------------------------------
                // Mterm46      = \sum(Bc) M[IJk](A,Bc) * T(Mk,Bc)
                // Mterm46      = \sum(Bc) M[IJk](Bc,A) * T(kM,Bc)
                // Mterm46      = \sum(Bc) T(kM,Bc)     * M[IJk](Bc,A) 
                SharedTensor2d Mtemp = std::make_shared<Tensor2d>("Mtemp", navirA * navirB, navirA);
                Mtemp->trans(Td_AAB);
                SharedTensor2d t2ABs = std::make_shared<Tensor2d>("t2ABs(iJ,Ab)", naoccB, naoccA, navirA, navirB);
                t2ABs->sort(2134, t2AB, 1.0, 0.0);
                // Mterm46(M,A) = \sum(Bc) T[k](M,Bc)   * M[IJk](Bc,A)
                SharedTensor2d Mterm46 = std::make_shared<Tensor2d>("Mterm46", naoccA, navirA);
                Mterm46->contract(false, false, naoccA, navirA, navirA * navirB, t2ABs, Mtemp,
                        (k * naoccA * navirA * navirB), 0, 1.0, 0.0);
                Mtemp.reset();
                t2ABs.reset();

                SharedTensor2d Mterm46tr = std::make_shared<Tensor2d>("Mterm46tr", navirA, naoccA);
                Mterm46tr->trans(Mterm46);
                Mterm46.reset();

                // M_OOVO(IJ,AM) += \sum(k) Mterm46tr(A,M)
                M_OOVO->axpy((size_t)(navirA * naoccA), 0, 1, Mterm46tr, (i * naoccA * navirA * naoccA) + (j * navirA * naoccA), 1, 1.0);
                Mterm46tr.reset();

                // UCCSD(T) M intr. Eq. (47)
                //------------------------------------------------------------------------------------------
                // Mterm47      = \sum(Bc) M[IJk](A,Bc) * T(Jm,Bc)
                // Mterm47      = \sum(Bc) M[IJk](Bc,A) * T(Jm,Bc)
                // Mterm47      = \sum(Bc) T(Jm,Bc)     * M[IJk](Bc,A) 
                Mtemp = std::make_shared<Tensor2d>("Mtemp", navirA * navirB, navirA);
                Mtemp->trans(Td_AAB);
                // Mterm47(m,A) = \sum(Bc) T[J](m,Bc)   * M[IJk](Bc,A)
                SharedTensor2d Mterm47 = std::make_shared<Tensor2d>("Mterm47(m,A)", naoccB, navirA);
                Mterm47->contract(false, false, naoccB, navirA, navirA * navirB, t2AB, Mtemp,
                        (j * naoccB * navirA * navirB), 0, 1.0, 0.0);
                Mtemp.reset();

                SharedTensor2d Mterm47tr = std::make_shared<Tensor2d>("Mterm47tr(A,m)", navirA, naoccB);
                Mterm47tr->trans(Mterm47);
                Mterm47.reset();

                // M_OoVo(Ik,Am) += \sum(J) Mterm47tr(A,m)
                M_OoVo->axpy((size_t)(navirA * naoccB), 0, 1, Mterm47tr, (i * naoccB * navirA * naoccB) + (k * navirA * naoccB), 1, 1.0);
                Mterm47tr.reset();

                // UCCSD(T) M intr. Eq. (50)
                //------------------------------------------------------------------------------------------
                // Mterm50       = \sum(B) Tc[IJk](A,Bc) * T(J,B)
                // Mterm50       = \sum(B) Tc[IJk](Bc,A) * T(J,B)
                // Mterm50       = \sum(B) T(J,B)        * Tc[IJk](Bc,A)
                SharedTensor2d Tc_AABtr = std::make_shared<Tensor2d>("Tc_AABtr", navirA * navirB, navirA);
                Tc_AABtr->trans(Tc_AAB);
                // Mterm50(cA) = \sum(B) T[J](B)        * Tc[IJk](B,cA)
                SharedTensor2d Mterm50 = std::make_shared<Tensor2d>("Mterm50(cA)", navirB, navirA);
                Mterm50->contract(false, false, 1, navirB * navirA, navirA, t1A, Tc_AABtr,
                        (j * navirA), 0, 1.0, 0.0);

                SharedTensor2d Mterm50t = std::make_shared<Tensor2d>("Mterm50t(A,c)", navirA, navirB);
                Mterm50t->trans(Mterm50);
                Mterm50.reset();

                // M_OoVv(IkAc) += \sum(J) Mterm50(Ac)
                SharedTensor2d M_OoVv = std::make_shared<Tensor2d>("M <Ij|Ab>", naoccA, naoccB, navirA, navirB);
                if (!first_write) M_OoVv->read(psio_, PSIF_DFOCC_DENS);
                M_OoVv->axpy((size_t)(navirA * navirB), 0, 1, Mterm50t, (i * naoccB * navirA * navirB) + (k * navirA * navirB), 1, 1.0);
                M_OoVv->write(psio_, PSIF_DFOCC_DENS);
                M_OoVv.reset();

                // UCCSD(T) M intr. Eq. (51)
                //------------------------------------------------------------------------------------------
                // Mterm51       = \sum(c) Tc[IJk](A,Bc) * T(k,c)
                // Mterm51(AB)   = \sum(c) Tc[IJk](AB,c)   * T[k](c)
                SharedTensor2d Mterm51 = std::make_shared<Tensor2d>("Mterm51(AB)", navirB, navirA);
                Mterm51->contract(false, false, navirA * navirA, 1, navirB, Tc_AAB, t1B,
                        0, (k * navirB), 1.0, 0.0);

                // M_OOVV(IJAB) += \sum(k) Mterm51(AB)
                SharedTensor2d M_OOVV = std::make_shared<Tensor2d>("M <IJ||AB>", naoccA, naoccA, navirA, navirA);
                M_OOVV->read_anti_symm(psio_, PSIF_DFOCC_DENS);
                M_OOVV->axpy((size_t)(navirA * navirA), 0, 1, Mterm51, (i * naoccA * navirA * navirA) + (j * navirA * navirA), 1, 1.0);
                M_OOVV->write_anti_symm(psio_, PSIF_DFOCC_DENS);
                M_OOVV.reset();

                Tc_AAB.reset();

                // UCCSD(T) Lambda Eq. (24) 1st term
                //------------------------------------------------------------------------------------------
                // read K[k](B,fE)
                K_i_AbC = std::make_shared<Tensor2d>("K[i] <A|bC>", navirA, navirB, navirA);
                K_i_AbC->myread(psio_, PSIF_DFOCC_IABC_BABA, (size_t)(k * navirB * navir2AA) * sizeof(double));

                SharedTensor2d lterm24_1 = std::make_shared<Tensor2d>("lterm24_1", navirA, navirA);
                SharedTensor2d lterm24_1tr = std::make_shared<Tensor2d>("lterm24_1", navirA, navirA);

                // lterm24_1(A,B) = \sum(fE) M[IJk](A,Ef) <kB|fE>
                // lterm24_1(A,B) = \sum(fE) M[IJk](A,Ef) K[k](B,fE)
                // lterm24_1(A,B) = \sum(fE) M[IJk](A,Ef) K[k](B,Ef)
                SharedTensor2d K_i_ACb = std::make_shared<Tensor2d>("K[i] <A|Cb>", navirA, navirA, navirB);
                K_i_ACb->sort3a(132, navirA, navirB, navirA, K_i_AbC, 1.0, 0.0);
                K_i_AbC.reset();
                // lterm24_1(A,B) = \sum(fE) M[IJk](A,Ef) K[k](Ef,B)
                lterm24_1->gemm(false, true, Td_AAB, K_i_ACb, 1.0, 0.0);
                K_i_ACb.reset();

                // tL2AA[IJ](AB) += P_(AB) \sum(k) lterm24_1(AB)
                tL2AA->axpy((size_t)navir2AA, 0, 1, lterm24_1, i * naoccA * navir2AA + j * navir2AA, 1, 1.0);
                lterm24_1tr->trans(lterm24_1);
                lterm24_1.reset();
                tL2AA->axpy((size_t)navir2AA, 0, 1, lterm24_1tr, i * naoccA * navir2AA + j * navir2AA, 1, -1.0);
                lterm24_1tr.reset();

                // UCCSD(T) Lambda Eq. (25)
                //------------------------------------------------------------------------------------------
                SharedTensor2d lterm25 = std::make_shared<Tensor2d>("lterm25", navir2AA, naoccA);
                // lterm25(AB,L) = -\sum(c) M[IJk](AB,c) K[Jk](c,L)
                lterm25->contract(false, false, navir2AA, naoccA, navirB, Td_AAB, K_IjaK,
                        0, (j * naoccB * naoccA * navirB) + (k * naoccA * navirB), -1.0, 0.0);

                SharedTensor2d lterm25_t = std::make_shared<Tensor2d>("lterm25_t", naoccA, navir2AA);
                lterm25_t->trans(lterm25);
                lterm25.reset();

                // tL2AA(IL,AB) += P_(IL) lterm25_t(L,AB)
                tL2AA->axpy((size_t)naoccA * navir2AA, 0, 1, lterm25_t, i * naoccA * navir2AA, 1, 1.0);

                for (int l = 0; l < naoccA; l++) {
                    long int li = ij_idxAA->get(l, i);
                    for (int a = 0; a < navirA; a++) {
                        for (int b = 0; b < navirA; b++) {
                            long int ab = ab_idxAA->get(a, b);
                            double val = lterm25_t->get(l, ab);
                            tL2AA->subtract(li, ab, val);
                        }
                    }
                }
                lterm25_t.reset();

                // UCCSD(T) Lambda Eq. (27) 1st term
                //------------------------------------------------------------------------------------------
                // read K[J](a,Bc)
                K_I_aBc = std::make_shared<Tensor2d>("K[I] <a|Bc>", navirB, navirA, navirB);
                K_I_aBc->myread(psio_, PSIF_DFOCC_IABC_ABAB, (size_t)(j * navirB * navirA * navirB) * sizeof(double));

                // lterm27(Ab) = \sum(Ef) M[IJk](A,Ef) <Jb|Ef>
                // lterm27(Ab) = \sum(Ef) M[IJk](A,Ef) K[J](b,Ef)
                SharedTensor2d lterm27 = std::make_shared<Tensor2d>("lterm27", navirA, navirB);
                lterm27->contract(false, true, navirA, navirB, navirA * navirB, Td_AAB, K_I_aBc, 1.0, 0.0);
                K_I_aBc.reset();

                // UCCSD(T) Lambda Eq. (27) 2nd term
                //------------------------------------------------------------------------------------------
                // read G[J](A,BC)
                G_I_ABC = std::make_shared<Tensor2d>("G[I] <A||BC>", navirA, navirA, navirA);
                G_I_ABC->myread(psio_, PSIF_DFOCC_IABC_AAAA, (size_t)(j * navirA * navir2AA) * sizeof(double));

                // lterm27(Ab) = -(1/2) \sum(EF) <JA||EF> M[IJk](EF,b)
                // lterm27(Ab) = -(1/2) \sum(EF) G[J](A,EF) M[IJk](EF,b)
                lterm27->contract(false, false, navirA, navirB, navir2AA, G_I_ABC, Td_AAB, -0.5, 1.0);
                G_I_ABC.reset();

                // tL2AB[Ik](Ab) += \sum(J) lterm27(Ab)
                tL2AB->axpy((size_t)navirA * navirB, 0, 1, lterm27, (i * naoccB * navirA * navirB) + (k * navirA * navirB), 1, 1.0);
                lterm27.reset();

                // UCCSD(T) Lambda Eq. (28)
                //------------------------------------------------------------------------------------------
                SharedTensor2d lterm28 = std::make_shared<Tensor2d>("lterm28 (l,cA)", naoccB, navirB, navirA);
                // lterm28(l,cA) = -\sum(B) <Jk|Bl> M[IJk](ABc)
                SharedTensor2d K_IjkA = std::make_shared<Tensor2d>("K_IjkA", naoccA, naoccB, naoccB, navirA);
                K_IjkA->sort(1243, K_IjAk, 1.0, 0.0);

                // lterm28(l,cA) = -\sum(B) <Jk|Bl> M[IJk](A,Bc)
                // sort
                // lterm28(l,cA) = -\sum(B) <Jk|lB> M[IJk](A,Bc)
                // lterm28(l,cA) = -\sum(B) [Jk](l,B) Mtemp[IJk](B,cA)
                Mtemp = std::make_shared<Tensor2d>("Mtemp <abC>", navirB, navirB, navirA);
                Mtemp->sort3a(231, navirA, navirA, navirB, Td_AAB, 1.0, 0.0);
                //lterm28->contract(false, true, naoccB, navirB * navirA, navirA, K_IjkA, Td_AAB,
                //        (j * naoccB * naoccB * navirA) + (k * naoccB * navirA), 0, -1.0, 0.0);
                lterm28->contract(false, false, naoccB, navirB * navirA, navirA, K_IjkA, Mtemp,
                        (j * naoccB * naoccB * navirA) + (k * naoccB * navirA), 0, -1.0, 0.0);
                K_IjkA.reset();

                SharedTensor2d lterm28_s = std::make_shared<Tensor2d>("lterm28_s (l,Ac)", naoccB, navirA, navirB);
                // (lcA) sort 132 --> (lAc)
                lterm28_s->sort3a(132, naoccB, navirB, navirA, lterm28, 1.0, 0.0);
                lterm28.reset();

                // tL2AB(Il,Ac) += \sum(Jk) lterm28(l,Ac)
                tL2AB->axpy((size_t)naoccB * navirA * navirB, 0, 1, lterm28_s, i * naoccB * navirA * navirB, 1, 1.0);
                lterm28_s.reset();

                // UCCSD(T) Lambda Eq. (29)
                //------------------------------------------------------------------------------------------
                SharedTensor2d lterm29 = std::make_shared<Tensor2d>("lterm29", naoccA, navirB * navirA);
                // lterm29(L,cA) = -0.5 \sum(B) M[IJk](ABc) <IJ||LB>
                SharedTensor2d G_IJKA = std::make_shared<Tensor2d>("G_IJKA", naoccA, naoccA, naoccA, navirA);
                G_IJKA->sort(1243, G_IJAK, -1.0, 0.0);
                // lterm29(L,cA) = -0.5 \sum(B) G_[IJ](L,B) M[IJk](A,Bc)
                // lterm29(L,cA) = -0.5 \sum(B) G_[IJ](L,B) Mtemp[IJk](B,cA)
                // already didt it above: Mtemp->sort3a(231, navirA, navirA, navirB, Td_AAB, 1.0, 0.0);
                //lterm29->contract(false, true, naoccA, navirB * navirA, navirA, G_IJKA, Td_AAB,
                //        (i * naoccA * naoccA * navirA) + (j * naoccA * navirA), 0, -0.5, 0.0);
                lterm29->contract(false, false, naoccA, navirB * navirA, navirA, G_IJKA, Mtemp,
                        (i * naoccA * naoccA * navirA) + (j * naoccA * navirA), 0, -0.5, 0.0);
                G_IJKA.reset();
                
                // tL2AB(Lk,Ac) += \sum(IJ) lterm29(L,cA)
                for (int l = 0; l < naoccA; l++) {
                    long int lk = ij_idxAB->get(l, k);
                    for (int a = 0; a < navirA; a++) {
                        for (int c = 0; c < navirB; c++) {
                            long int ac = ab_idxAB->get(a, c);
                            long int ca = ab_idxBA->get(c, a);
                            double val = lterm29->get(l, ca);
                            tL2AB->add(lk, ac, val);
                        }
                    }
                }
                lterm29.reset();

                // UCCSD(T) M intr. Eq. (48)
                //------------------------------------------------------------------------------------------
                // Mterm48       = \sum(c) M[IJk](A,Bc) * T(Jk,Dc)
                // Mterm48       = \sum(c) M[IJk](A,Bc) * T(kJ,cD)
                // Mterm48(AB,D) = \sum(c) M[IJk](AB,c) * T[kJ](c,D)
                SharedTensor2d Mterm48 = std::make_shared<Tensor2d>("Mterm48(AB,D)", navir2AA, navirA);
                Mterm48->myread(psio_, PSIF_DFOCC_MIABC_AAAA, (size_t)(i * navirA * navir2AA) * sizeof(double));
                Mterm48->contract(false, false, navir2AA, navirA, navirB, Td_AAB, t2BA,
                        0, (k * naoccA * navirB * navirA) + (j * navirB * navirA), 1.0, 1.0);

                Mterm48->mywrite(psio_, PSIF_DFOCC_MIABC_AAAA, (size_t)(i * navirA * navir2AA) * sizeof(double));
                Mterm48.reset();

                // UCCSD(T) M intr. Eq. (49)
                //------------------------------------------------------------------------------------------
                // Mterm49       = \sum(B) M[IJk](A,Bc) * T(Jk,Bd)
                // Mterm49       = \sum(B) M[IJk](Bc,A) * T(kJ,dB)
                Mtemp = std::make_shared<Tensor2d>("Mtemp", navirA * navirB, navirA);
                Mtemp->trans(Td_AAB);
                // Mterm49(d,cA) = \sum(B) T[kJ](d,B)    * M[IJk](B,cA)
                SharedTensor2d Mterm49 = std::make_shared<Tensor2d>("Mterm49(d,cA)", navirB, navirB, navirA);
                Mterm49->contract(false, false, navirB, navirB * navirA, navirA, t2BA, Mtemp,
                        (k * naoccA * navirB * navirA) + (j * navirB * navirA), 0, 1.0, 0.0);
                Mtemp.reset();

                Mterm49s->sort3a(321, navirB, navirB, navirA, Mterm49, 1.0, 1.0);
                Mterm49.reset();

                // UCCSD(T) M intr. Eq. (45)
                //------------------------------------------------------------------------------------------
                // Mterm45       = (1/2) *  \sum(A) M[IJk](A,Bc) * T(IJ,AD)
                // Mterm45       = (1/2) * -\sum(A) M[IJk](A,Bc) * T(IJ,DA)
                // Mterm45(D,Bc) = (1/2) * -\sum(A) T[IJ](D,A)     * M[IJk](A,Bc)
                SharedTensor2d Mterm45 = std::make_shared<Tensor2d>("Mterm45", navirA, navirA, navirB);
                Mterm45->contract(false, false, navirA, navirA * navirB, navirA, t2AA, Td_AAB,
                        (i * naoccA * navir2AA) + (j * navir2AA), 0, -0.5, 0.0);

                SharedTensor2d Mterm45s = std::make_shared<Tensor2d>("Mterm45tr(c,BD)", navirB, navirA, navirA);
                Mterm45s->myread(psio_, PSIF_DFOCC_MIABC_BBAA, (size_t)(k * navirB * navir2AA) * sizeof(double));
                Mterm45s->sort3a(321, navirA, navirA, navirB, Mterm45, 1.0, 1.0);
                Mterm45.reset();
                Mterm45s->mywrite(psio_, PSIF_DFOCC_MIABC_BBAA, (size_t)(k * navirB * navir2AA) * sizeof(double));
                Mterm45s.reset();

                Td_AAB.reset();

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
        Mterm49s->mywrite(psio_, PSIF_DFOCC_MIABC_AABB, (bool)i);
        Mterm49s->zero();
    } // i

    Mterm49s.reset();

    outfile->Printf("\tAAB (T) energy                     : % 20.14f\n", sumAAB);
    E_t += sumAAB;

    // reset all AAB things
    G_IJAK.reset();
    G_IJAB.reset();
    t2AA.reset();
    K_I_aBc.reset();
    K_IjAb.reset();
    timer_off("CCSD(T)-HM-AAB");

    timer_on("CCSD(T)-HM-ABB");
    //==================================================
    //======================= ABB ======================
    //==================================================
    // N10 : beginning of ABB

    // malloc W_ABB[Ijk](Abc)
    SharedTensor2d W_ABB = std::make_shared<Tensor2d>("X[Ijk](A,cb)", navirA, navirB, navirB);

    // read T2BB
    t2BB = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    t2BB->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

    // progress counter
    start = std::time(nullptr);
    stop = std::time(nullptr);
    ind = 0;
    next_print = step_print;

    Nijk = naoccA * naoccB * naoccB;
    outfile->Printf("\n\tNumber of ijk combinations for ABB: %i \n", Nijk);

    // N11 : main loop of ABB
    double sumABB = 0.0;
    for (long int i = 0; i < naoccA; i++) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);

        for(long int j = 0; j < naoccB; j++) {
            long int ij = ij_idxAB->get(i, j);
            double Dij = Di + FockB->get(j + nfrzc, j + nfrzc);

            for (long int k = 0; k < naoccB; k++) {
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

                double Dijk = Dij + FockB->get(k + nfrzc, k + nfrzc);

                // read <ij||ab>
                G_ijab = std::make_shared<Tensor2d>("G <ij||ab>", naoccB, naoccB, navirB, navirB);
                G_ijab->read(psio_, PSIF_DFOCC_IJAB_BBBB);

                K_IjAb = std::make_shared<Tensor2d>("K <Ij|Ab>", naoccA, naoccB, navirA, navirB);
                K_IjAb->read(psio_, PSIF_DFOCC_IJAB_ABAB);

                // malloc Tc[Ijk](Abc)
                SharedTensor2d Tc_ABB = std::make_shared<Tensor2d>("Tc[Ijk](A,bc)", navirA, navirB, navirB);
                // malloc Td[Ijk](Abc)
                SharedTensor2d Td_ABB = std::make_shared<Tensor2d>("Td[Ijk](A,bc)", navirA, navirB, navirB);

                double Wijkabc, Vijkabc;
//#pragma omp parallel for private(Wijkabc, Vijkabc) reduction(+ : sumABB)
                for (long int a = 0; a < navirA; a++) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    double LIA = 0;
                    for (long int b = 0; b < navirB; b++) {
                        double Dijkab = Dijka - FockB->get(b + noccB, b + noccB);
                        long int ab = ab_idxAB->get(a, b);
                        long int ba = ab_idxBA->get(b, a);
                        for (long int c = 0; c < navirB; c++) {
                            long int ac = ab_idxAB->get(a, c);
                            long int ca = ab_idxBA->get(c, a);
                            long int cb = ab_idxBB->get(c, b);
                            long int bc = ab_idxBB->get(b, c);

                            double Dijkabc = Dijkab - FockB->get(c + noccB, c + noccB);

                            // W_ABB[Ijk](Abc)
                            Wijkabc = W_ABB->get(a, cb);

                            // V_ABB = t_Ij^Ab * f_kc + t_k^c * <Ij|Ab>
                            //       - t_Ik^Ab * f_jc - t_j^c * <Ik|Ab>
                            //       - t_Ij^Ac * f_kb - t_k^b * <Ij|Ac>
                            //       + t_Ik^Ac * f_jb + t_j^b * <Ik|Ac>
                            //       + t_kj^cb * f_I^A + t_I^A * <kj||cb>
                            Vijkabc = t2AB->get(ij, ab) * FockB->get(nfrzc + k, noccB + c) + t1B->get(k, c) * K_IjAb->get(ij, ab)
                                    - t2AB->get(ik, ab) * FockB->get(nfrzc + j, noccB + c) - t1B->get(j, c) * K_IjAb->get(ik, ab)
                                    - t2AB->get(ij, ac) * FockB->get(nfrzc + k, noccB + b) - t1B->get(k, b) * K_IjAb->get(ij, ac)
                                    + t2AB->get(ik, ac) * FockB->get(nfrzc + j, noccB + b) + t1B->get(j, b) * K_IjAb->get(ik, ac)
                                    + t2BB->get(kj, cb) * FockA->get(nfrzc + i, noccA + a) + t1A->get(i, a) * G_ijab->get(kj, cb);

                            double tc_ijkabc = Wijkabc / Dijkabc;
                            Tc_ABB->set(a, bc, tc_ijkabc);

                            LIA += (0.25) * G_ijab->get(jk, bc) * tc_ijkabc;

                            double Lia = K_IjAb->get(ij, ab) * tc_ijkabc;
                            tL1B->add(k, c, Lia);

                            double td_ijkabc = Vijkabc / Dijkabc;
                            Td_ABB->set(a, bc, td_ijkabc);

                            // N12 : energy of ABB
                            sumABB += 0.25 * (Wijkabc + Vijkabc) * Wijkabc / Dijkabc;

                        } // c
                    } // b
                    tL1A->add(i, a, LIA);
                } // a

                // tijkabc = t(c)ijkabc + t(d)ijkabc
                Td_ABB->axpy(Tc_ABB, 1.0);

                // UCCSD(T) PDM Eq. (3), (4), (7) and (10)
                for (int a = 0; a < navirA; a++) {
                    double G1c_ABB_contrib = 0;
                    for (int b = 0; b < navirB; b++) {
                        double G1c_ABB_contrib3 = 0;
                        for (int c = 0; c < navirB; c++) {
                            long int bc = ab_idxBB->get(b, c);
                            G1c_ABB_contrib += (1.0/4.0) * Td_ABB->get(a, bc) * Tc_ABB->get(a, bc);

                            double G1c_ABB_contrib2 = (1.0/4.0) * Td_ABB->get(a, bc) * Tc_ABB->get(a, bc);
                            G1c_iiB->subtract(k, G1c_ABB_contrib2);
                            G1c_aaB->add(c, G1c_ABB_contrib2);

                            G1c_ABB_contrib3 += (1.0/4.0) * Td_ABB->get(a, bc) * Tc_ABB->get(a, bc);
                        }
                        G1c_iiB->subtract(j, G1c_ABB_contrib3);
                        G1c_aaB->add(b, G1c_ABB_contrib3);
                    }
                    G1c_iiA->subtract(i, G1c_ABB_contrib);
                    G1c_aaA->add(a, G1c_ABB_contrib);
                }

                // Mijkabc = tijkabc + t(c)ijkabc
                Td_ABB->axpy(Tc_ABB, 1.0);

                // UCCSD(T) M intr. Eq. (52)
                //------------------------------------------------------------------------------------------
                // Mterm52      = \sum(Ab) M[Ijk](A,bc) * T(Im,Ab)
                // Mterm52      = \sum(Ab) T(Im,Ab)     * M[Ijk](A,bc) 
                // Mterm52(m,c) = \sum(Ab) T[I](m,Ab)   * M[Ijk](Ab,c) 
                SharedTensor2d Mterm52 = std::make_shared<Tensor2d>("Mterm52", naoccB, navirB);
                Mterm52->contract(false, false, naoccB, navirB, navirA * navirB, t2AB, Td_ABB,
                        (i * naoccB * navirA * navirB), 0, 1.0, 0.0);

                SharedTensor2d Mterm52tr = std::make_shared<Tensor2d>("Mterm52tr(c,m)", navirB, naoccB);
                Mterm52tr->trans(Mterm52);
                Mterm52.reset();

                // M_oovo(kj,cm) += \sum(I) Mterm52(c,m)
                M_oovo->axpy((size_t)(naoccB * navirB), 0, 1, Mterm52tr, (k * naoccB * navirB * naoccB) + (j * navirB * naoccB), 1, 1.0);
                Mterm52tr.reset();

                // UCCSD(T) M intr. Eq. (53)
                // Crawford i
                //------------------------------------------------------------------------------------------
                // Mterm53      =  \sum(Ac) M[Ijk](A,bc) * T(Mk,Ac)
                // Mterm53      = -\sum(Ac) M[Ijk](A,cb) * T(kM,Ac)
                SharedTensor2d t2ABs = std::make_shared<Tensor2d>("t2ABs(iJ,Ab)", naoccB, naoccA, navirA, navirB);
                t2ABs->sort(2134, t2AB, 1.0, 0.0);
                // Mterm53      = -\sum(Ac) T(kM,Ac)     * M[Ijk](A,cb)
                // Mterm53(M,b) = -\sum(Ac) T[k](M,Ac)   * M[Ijk](Ac,b)
                SharedTensor2d Mterm53 = std::make_shared<Tensor2d>("Mterm53", naoccA, navirB);
                Mterm53->contract(false, false, naoccA, navirB, navirA * navirB, t2ABs, Td_ABB,
                        (k * naoccA * navirA * navirB), 0, -1.0, 0.0);
                t2ABs.reset();

                SharedTensor2d Mterm53tr = std::make_shared<Tensor2d>("Mterm53tr(b,M)", navirB, naoccA);
                Mterm53tr->trans(Mterm53);
                Mterm53.reset();

                // M_oOvO[jI](bM) += \sum(k) Mterm53tr(bM)
                M_oOvO->axpy((size_t)(navirB * naoccA), 0, 1, Mterm53tr, (j * naoccA * navirB * naoccA) + (i * navirB * naoccA), 1, 1.0);
                Mterm53tr.reset();

                // UCCSD(T) M intr. Eq. (55)
                //------------------------------------------------------------------------------------------
                // Mterm55       = \sum(c) M[Ijk](A,bc) * T(Ik,Dc)
                // Mterm55       = \sum(c) M[Ijk](A,bc) * T(kI,cD)
                // Mterm55(AbD)  = \sum(c) M[Ijk](Ab,c) * T[kI](c,D)
                SharedTensor2d Mterm55 = std::make_shared<Tensor2d>("Mterm55(A,bD)", navirA, navirB, navirA);
                Mterm55->contract(false, false, navirA * navirB, navirA, navirB, Td_ABB, t2BA,
                        0, (k * naoccA * navirB * navirA) + (i * navirB * navirA), 1.0, 0.0);

                SharedTensor2d Mterm55s = std::make_shared<Tensor2d>("Mterm55tr(b,AD)", navirB, navirA, navirA);
                Mterm55s->myread(psio_, PSIF_DFOCC_MIABC_BBAA, (size_t)(j * navirB * navir2AA) * sizeof(double));
                Mterm55s->sort3a(213, navirA, navirB, navirA, Mterm55, 1.0, 1.0);
                Mterm55.reset();
                Mterm55s->mywrite(psio_, PSIF_DFOCC_MIABC_BBAA, (size_t)(j * navirB * navir2AA) * sizeof(double));
                Mterm55s.reset();

                // UCCSD(T) M intr. Eq. (56)
                //------------------------------------------------------------------------------------------
                // Mterm56        = (0.5) * -\sum(bc) M[Ijk](A,bc) * T(km,bc)
                // Mterm56        = (0.5) * -\sum(bc) T(km,bc)     * M[Ijk](bc,A)
                // Mterm56(m,A)   = (0.5) * -\sum(bc) T[k](m,bc)     * M[Ijk](bc,A)
                SharedTensor2d Mterm56 = std::make_shared<Tensor2d>("Mterm56(m,A)", naoccB, navirA);
                Mterm56->contract(false, true, naoccB, navirA, navir2BB, t2BB, Td_ABB,
                        (k * naoccB * navir2BB), 0, -0.5, 0.0);

                SharedTensor2d Mterm56tr = std::make_shared<Tensor2d>("Mterm56tr(A,m)", navirA, naoccB);
                Mterm56tr->trans(Mterm56);
                Mterm56.reset();

                // M_OoVo(Ij,Am) += \sum(k) Mterm56tr(A,m)
                M_OoVo->axpy((size_t)(navirA * naoccB), 0, 1, Mterm56tr, (i * naoccB * navirA * naoccB) + (j * navirA * naoccB), 1, 1.0);
                Mterm56tr.reset();

                // UCCSD(T) M intr. Eq. (57)
                //------------------------------------------------------------------------------------------
                // Mterm57       =  \sum(c) M[Ijk](A,bc) * T(jk,dc)
                // Mterm57       = -\sum(c) M[Ijk](A,bc) * T(jk,cd)
                // Mterm57(Ab,d) = -\sum(c) M[Ijk](Ab,c) * T[jk](c,d)
                SharedTensor2d Mterm57 = std::make_shared<Tensor2d>("Mterm57(Ab,d)", navirA, navirB, navirB);
                Mterm57->myread(psio_, PSIF_DFOCC_MIABC_AABB,(size_t)(i * navirA * navir2BB) * sizeof(double));
                Mterm57->contract(false, false, navirA * navirB, navirB, navirB, Td_ABB, t2BB,
                        0, (j * naoccB * navirB * navirB) + (k * navirB * navirB), -0.5, 1.0);
                Mterm57->mywrite(psio_, PSIF_DFOCC_MIABC_AABB,(size_t)(i * navirA * navir2BB) * sizeof(double));
                Mterm57.reset();

                // UCCSD(T) M intr. Eq. (58)
                //------------------------------------------------------------------------------------------
                // Mterm58(bc) = \sum(A) T[I](A) * Tc[Ijk](A,bc)
                SharedTensor2d Mterm58 = std::make_shared<Tensor2d>("Mterm58", navirB, navirB);
                Mterm58->contract(false, false, 1, navir2BB, navirA, t1A, Tc_ABB,
                        (i * navirA), 0, 1.0, 0.0);

                SharedTensor2d Mterm58tr = std::make_shared<Tensor2d>("Mterm58tr(cb)", navirB, navirB);
                Mterm58tr->trans(Mterm58);
                Mterm58.reset();

                // M_oovv(kj,cb) += \sum(I) Mterm58(bc)
                SharedTensor2d M_oovv = std::make_shared<Tensor2d>("M <ij||ab>", naoccB, naoccB, navirB, navirB);
                M_oovv->read_anti_symm(psio_, PSIF_DFOCC_DENS);
                M_oovv->axpy((size_t)(navirB * navirB), 0, 1, Mterm58tr, (k * naoccB * navirB * navirB) + (j * navirB * navirB), 1, 1.0);
                M_oovv->write_anti_symm(psio_, PSIF_DFOCC_DENS);
                M_oovv.reset();
                Mterm58tr.reset();

                // UCCSD(T) M intr. Eq. (59)
                //------------------------------------------------------------------------------------------
                // Mterm59       = \sum(c) Tc[Ijk](A,bc) * T(k,c)
                // Mterm59(Ab)   = \sum(c) Tc[Ijk](Ab,c) * T[k](c) 
                SharedTensor1d Mterm59 = std::make_shared<Tensor1d>("Mterm59(Ab)", navirA * navirB);
                Mterm59->gemv(false, navirA * navirB, navirB, Tc_ABB, t1B,
                        0, (k * navirB), 1.0, 0.0);

                // M_OoVv(IjAb) += \sum(k) Mterm59(Ab)
                SharedTensor2d M_OoVv = std::make_shared<Tensor2d>("M <Ij|Ab>", naoccA, naoccB, navirA, navirB);
                M_OoVv->read(psio_, PSIF_DFOCC_DENS);
                M_OoVv->axpy((size_t)(navirA * navirB), 0, 1, Mterm59, (i * naoccB * navirA * navirB) + (j * navirA * navirB), 1, 1.0);
                M_OoVv->write(psio_, PSIF_DFOCC_DENS);
                M_OoVv.reset();
                Mterm59.reset();

                // UCCSD(T) Lambda Eq. (31)
                //------------------------------------------------------------------------------------------
                // read K[I](a,Bc)
                K_I_aBc = std::make_shared<Tensor2d>("K[I] <a|Bc>", navirB, navirA, navirB);
                K_I_aBc->myread(psio_, PSIF_DFOCC_IABC_ABAB, (size_t)(i * navirB * navirA * navirB) * sizeof(double));

                SharedTensor2d lterm31 = std::make_shared<Tensor2d>("lterm31", navirB, navirB);
                SharedTensor2d lterm31_tr = std::make_shared<Tensor2d>("lterm31", navirB, navirB);
                // lterm31(ba) = \sum(Fe) <Ib|Fe> M[Ijk](Fe,a)
                // lterm31(ba) = \sum(Fe) K[I](b,Fe) M[Ijk](Fe,a)
                lterm31->contract(false, false, navirB, navirB, navirA * navirB, K_I_aBc, Td_ABB, 1.0, 0.0);
                K_I_aBc.reset();

                // tL2BB[kj](ab) += P_(ab) \sum(I) lterm31(ab)
                tL2BB->axpy((size_t)navirB * navirB, 0, 1, lterm31, (k * naoccB * navir2BB) + (j * navir2BB), 1, -1.0);
                lterm31_tr->trans(lterm31);
                lterm31.reset();
                tL2BB->axpy((size_t)navirB * navirB, 0, 1, lterm31_tr, (k * naoccB * navir2BB) + (j * navir2BB), 1, 1.0);
                lterm31_tr.reset();

                // UCCSD(T) Lambda Eq. (32)
                //------------------------------------------------------------------------------------------
                SharedTensor2d lterm32 = std::make_shared<Tensor2d>("lterm32", navirB * navirB, naoccB);
                // lterm32(ab,l) = \sum(C) M[Ijk](Cab) K[Ik](C,l)
                // lterm32(ab,l) = \sum(C) Mtemp[Ijk](ab,C) K[Ik](C,l)
                SharedTensor2d Mtemp = std::make_shared<Tensor2d>("Mtemp <abC>", navirB, navirB, navirA);
                Mtemp->sort3a(231, navirA, navirB, navirB, Td_ABB, 1.0, 0.0);
                //lterm32->contract(true, false, navir2BB, naoccB, navirA, Td_ABB, K_IjAk,
                //        0, (i * naoccB * navirA * naoccB) + (k * navirA * naoccB), -1.0, 0.0);
                lterm32->contract(false, false, navir2BB, naoccB, navirA, Mtemp, K_IjAk,
                        0, (i * naoccB * navirA * naoccB) + (k * navirA * naoccB), -1.0, 0.0);
                Mtemp.reset();

                SharedTensor2d lterm32_t = std::make_shared<Tensor2d>("lterm32_t", naoccB, navir2BB);
                lterm32_t->trans(lterm32);
                lterm32.reset();

                // tL2BB[j](l,ab) += P_(jl) \sum(Ik) lterm32(l,ab)
                tL2BB->axpy((size_t)naoccB * navir2BB, 0, 1, lterm32_t, j * naoccB * navir2BB, 1, 1.0);

                for (int l = 0; l < naoccB; l++) {
                    long int lj = ij_idxBB->get(l, j);
                    for (int a = 0; a < navirB; a++) {
                        for (int b = 0; b < navirB; b++) {
                            long int ab = ab_idxBB->get(a, b);
                            double val = lterm32_t->get(l, ab);
                            tL2BB->subtract(lj, ab, val);
                        }
                    }
                }
                lterm32_t.reset();

                // UCCSD(T) Lambda Eq. (33) 1st term
                //------------------------------------------------------------------------------------------
                // read G[k](a,bc)
                G_i_abc = std::make_shared<Tensor2d>("G[i] <a||bc>", navirB, navirB, navirB);
                G_i_abc->myread(psio_, PSIF_DFOCC_IABC_BBBB, (size_t)(k * navirB * navir2BB) * sizeof(double));

                SharedTensor2d lterm33_1 = std::make_shared<Tensor2d>("lterm33_1", navirA, navirB);
                // lterm33_1(A,b) = 0.5 \sum(ef) M[Ijk](A,ef] <bk||ef> 
                // lterm33_1(A,b) = -0.5 \sum(ef) M[Ijk](A,ef] <kb||ef> 
                // lterm33_1(A,b) = -0.5 \sum(ef) M[Ijk](A,ef) G_[k](ef,b)
                lterm33_1->contract(false, true, navirA, navirB, navir2BB, Td_ABB, G_i_abc, -0.5, 0.0);
                G_i_abc.reset();

                // tL2AB[Ij](Ab) += \sum(k) lterm33_1(Ab)
                tL2AB->axpy((size_t)navirA * navirB, 0, 1, lterm33_1, (i * naoccB * navirA * navirB) + (j * navirA * navirB), 1, 1.0);
                lterm33_1.reset();

                // UCCSD(T) Lambda Eq. (33) 2nd term
                //------------------------------------------------------------------------------------------
                // read K[k](A,bC)
                K_i_AbC = std::make_shared<Tensor2d>("K[i] <A|bC>", navirA, navirB, navirA);
                K_i_AbC->myread(psio_, PSIF_DFOCC_IABC_BABA, (size_t)(k * navirA * navirB * navirA) * sizeof(double));

                SharedTensor2d K_i_ACb = std::make_shared<Tensor2d>("K_i_ACb", navirA, navirA, navirB);
                K_i_ACb->sort3a(132, navirA, navirB, navirA, K_i_AbC, 1.0, 0.0);
                K_i_AbC.reset();

                SharedTensor2d lterm33_2 = std::make_shared<Tensor2d>("lterm33_2", navirA, navirB);
                // lterm33_2(A,b) = -\sum(Ef) [k](A,Ef) M[Ijk](Ef,b)
                lterm33_2->contract(false, false, navirA, navirB, navirA * navirB, K_i_ACb, Td_ABB, -1.0, 0.0);
                K_i_ACb.reset();

                // tL2AB[Ij](Ab) += \sum(k) lterm33_2(Ab)
                tL2AB->axpy((size_t)navirA * navirB, 0, 1, lterm33_2, (i * naoccB * navirA * navirB) + (j * navirA * navirB), 1, 1.0);
                lterm33_2.reset();

                // UCCSD(T) Lambda Eq. (34)
                //------------------------------------------------------------------------------------------
                SharedTensor2d lterm34 = std::make_shared<Tensor2d>("lterm34", navirA * navirB, naoccB);
                // lterm34(Ab,l) = \sum(c) M[Ijk](Ab,c) <jk||cl>
                // lterm34(Ab,l) = \sum(c) M[Ijk](Ab,c) G[jk](c,l)
                lterm34->contract(false, false, navirA * navirB, naoccB, navirB, Td_ABB, G_ijak,
                        0, j * naoccB * navirB * naoccB + k * navirB * naoccB, 0.5, 0.0);

                // (Ab,l) -> (l,Ab)
                SharedTensor2d lterm34_t = std::make_shared<Tensor2d>("lterm34_t", naoccB, navirA, navirB);
                lterm34_t->trans(lterm34);
                lterm34.reset();

                // tL2AB[I](l,Ab) += \sum(jk) lterm34_t(l,Ab)
                tL2AB->axpy((size_t)naoccB * navirA * navirB, 0, 1, lterm34_t, i * naoccB * navirA * navirB, 1, 1.0);
                lterm34_t.reset();

                // UCCSD(T) Lambda Eq. (35)
                //------------------------------------------------------------------------------------------
                SharedTensor2d lterm35 = std::make_shared<Tensor2d>("lterm35", navirA * navirB, naoccA);

                // lterm35(Ab,L) = -\sum(c) M[Ijk](Ab,c) <Ik|Lc>
                // lterm35(Ab,L) = -\sum(c) M[Ijk](Ab,c) [Ik](c,L)
                lterm35->contract(false, false, navirA * navirB, naoccA, navirB, Td_ABB, K_IjaK,
                        0, (i * naoccB * navirB * naoccA) + (k * navirB * naoccA), -1.0, 0.0);

                // tL2AB(Lj,Ab) += \sum(Ik) lterm35(Ab,L)
                for (int l = 0; l < naoccA; l++) {
                    long int lj = ij_idxAB->get(l, j);
                    for (int a = 0; a < navirA; a++) {
                        for (int b = 0; b < navirB; b++) {
                            long int ab = ab_idxAB->get(a, b);
                            double val = lterm35->get(ab, l);
                            tL2AB->add(lj, ab, val);
                        }
                    }
                }
                lterm35.reset();

                // UCCSD(T) M intr. Eq. (54)
                //------------------------------------------------------------------------------------------
                // Mterm54       = -\sum(A) M[Ijk](Ab,c) * T(Ik,Ad)
                // Mterm54       = -\sum(A) T(kI,dA)     * M[Ijk](Ab,c)
                // Mterm54(d,bc) = -\sum(A) T[kI](d,A)   * M[Ijk](A,bc)
                SharedTensor2d Mterm54 = std::make_shared<Tensor2d>("Mterm54", navirB, navirB, navirB);
                Mterm54->contract(false, false, navirB, navir2BB, navirA, t2BA, Td_ABB,
                        (k * naoccA * navirB * navirA) + (i * navirB * navirA), 0, -1.0, 0.0);

                SharedTensor2d Mterm54s = std::make_shared<Tensor2d>("Mterm54tr(c,bd)", navirB, navirB, navirB);
                Mterm54s->myread(psio_, PSIF_DFOCC_MIABC_BBBB, (size_t)(j * navirB * navir2BB) * sizeof(double));
                Mterm54s->sort3a(321, navirB, navirB, navirB, Mterm54, 1.0, 1.0);
                Mterm54.reset();
                Mterm54s->mywrite(psio_, PSIF_DFOCC_MIABC_BBBB, (size_t)(j * navirB * navir2BB) * sizeof(double));
                Mterm54s.reset();

                Tc_ABB.reset();
                Td_ABB.reset();

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

    outfile->Printf("\tABB (T) energy                     : % 20.14f\n\n", sumABB);
    E_t += sumABB;

    // set energy
    Eccsd_t = Eccsd + E_t;

    // reset all ABB things
    W_ABB.reset();
    G_ijak.reset();
    K_IjAk.reset();
    K_IjKa.reset();
    K_IjaK.reset();
    G_ijab.reset();
    t2AB.reset();
    t2BA.reset();
    t2BB.reset();

    timer_off("CCSD(T)-HM-ABB");

    // remove files
    remove_binary_file(PSIF_DFOCC_IABC_AAAA);
    remove_binary_file(PSIF_DFOCC_IABC_BBBB);
    remove_binary_file(PSIF_DFOCC_IABC_BABA);
    remove_binary_file(PSIF_DFOCC_IABC_ABAB);

    tL1A->write(psio_, PSIF_DFOCC_AMPS);
    tL1A.reset();
    tL1B->write(psio_, PSIF_DFOCC_AMPS);
    tL1B.reset();

    tL2AA->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    tL2AA.reset();
    tL2BB->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
    tL2BB.reset();
    tL2AB->write(psio_, PSIF_DFOCC_AMPS);
    tL2AB.reset();

    M_OOVO->write(psio_, PSIF_DFOCC_DENS);
    M_oovo->write(psio_, PSIF_DFOCC_DENS);
    M_OoVo->write(psio_, PSIF_DFOCC_DENS);
    M_oOvO->write(psio_, PSIF_DFOCC_DENS);

} // uccsd_triples_hm()

}  // namespace dfoccwave
}  // namespace psi
