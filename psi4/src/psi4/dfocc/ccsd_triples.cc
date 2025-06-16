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

#include <ctime>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/libqt/qt.h"

#include "defines.h"
#include "dfocc.h"

namespace psi {
namespace dfoccwave {

void DFOCC::ccsd_canonic_triples() {
    // defs
    SharedTensor2d K, L, M, I, J, T, U, Tau, W, X, Y, Z;
    SharedTensor2d V, J1, J2, J3, Jt;
    SharedTensor1d Eijk;
    long int Nijk;

    // progress counter
    std::time_t stop, start = std::time(nullptr);
    long int ind = 0;
    double step_print = 10.0;
    double next_print = step_print;

    // Find number of unique ijk combinations (i>=j>=k)
    /*
    Nijk = 0;
    for(long int i = 0 ; i < naoccA; ++i){
        for(long int j = 0 ; j <= i; ++j){
            for(long int k = 0 ; k <= j; ++k){
                Nijk++;
            }
        }
    }
    */
    Nijk = naoccA * (naoccA + 1) * (naoccA + 2) / 6;
    outfile->Printf("\tNumber of ijk combinations: %i \n", Nijk);

    // Malloc Eijk
    // Eijk = std::make_shared<Tensor1d>("Eijk", Nijk);

    // Memory: 2*O^2V^2 + 5*V^3 + O^3V + V^2N + V^3/2

    // Read t2 amps
    t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->sort(1324, t2, 1.0, 0.0);
    t2.reset();

    // Form (ij|ka)
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    M->read(psio_, PSIF_DFOCC_INTS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, K, M, 1.0, 0.0);
    K.reset();
    I = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Form (ia|jb)
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, M, M, 1.0, 0.0);

    // B(iaQ)
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ);
    L = M->transpose();
    M.reset();

    // malloc W[ijk](abc)
    W = std::make_shared<Tensor2d>("W[IJK] <AB|C>", navirA * navirA, navirA);
    V = std::make_shared<Tensor2d>("V[IJK] <BA|C>", navirA * navirA, navirA);
    J1 = std::make_shared<Tensor2d>("J[I] <AB|E>", navirA * navirA, navirA);
    J2 = std::make_shared<Tensor2d>("J[J] <AB|E>", navirA * navirA, navirA);
    J3 = std::make_shared<Tensor2d>("J[K] <AB|E>", navirA * navirA, navirA);

    // B(Q,ab)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    K->read(psio_, PSIF_DFOCC_INTS);
    Jt = std::make_shared<Tensor2d>("J[I] <A|B>=C", navirA, ntri_abAA);

    // main loop
    E_t = 0.0;
    double sum = 0.0;
    for (long int i = 0; i < naoccA; ++i) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);

        // Compute J[i](a,bc) = (ia|bc) = \sum(Q) B[i](aQ) * B(Q,bc)
        Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, i * navirA * nQ, 0, 1.0, 0.0);
        J1->expand23(navirA, navirA, navirA, Jt);

        for (long int j = 0; j <= i; ++j) {
            long int ij = ij_idxAA->get(i, j);
            double Dij = Di + FockA->get(j + nfrzc, j + nfrzc);

            // Compute J[j](a,bc) = (ja|bc) = \sum(Q) B[j](aQ) * B(Q,bc)
            Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, j * navirA * nQ, 0, 1.0, 0.0);
            J2->expand23(navirA, navirA, navirA, Jt);

            for (long int k = 0; k <= j; ++k) {
                long int ik = ij_idxAA->get(i, k);
                long int jk = ij_idxAA->get(j, k);
                // Compute J[k](a,bc) = (ka|bc) = \sum(Q) B[k](aQ) * B(Q,bc)
                Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, k * navirA * nQ, 0, 1.0, 0.0);
                J3->expand23(navirA, navirA, navirA, Jt);

                // W[ijk](ab,c) = \sum(e) t_jk^ec (ia|be) (1+)
                // W[ijk](ab,c) = \sum(e) J[i](ab,e) T[jk](ec)
                W->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (j * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ab,c) -= \sum(m) t_im^ab <jk|mc> (1-)
                // W[ijk](ab,c) -= \sum(m) T[i](m,ab) I[jk](mc)
                W->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);

                // W[ijk](ac,b) = \sum(e) t_kj^eb (ia|ce) (2+)
                // W[ijk](ac,b) = \sum(e) J[i](ac,e) T[kj](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (k * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ac,b) -= \sum(m) t_im^ac <kj|mb> (2-)
                // W[ijk](ac,b) -= \sum(m) T[i](m,ac) I[kj](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA * navirA + b, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](ba,c) = \sum(e) t_ik^ec (jb|ae) (3+)
                // W[ijk](ba,c) = \sum(e) J[j](ba,e) T[ik](ec)
                V->contract(false, false, navirA * navirA, navirA, navirA, J2, T, 0,
                            (i * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ba,c) -= \sum(m) t_jm^ba <ik|mc> (3-)
                // W[ijk](ba,c) -= \sum(m) T[j](m,ba) I[ik](mc)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a * navirA, 1, V,
                                a * navirA * navirA + b * navirA, 1, 1.0);
                    }
                }

                // W[ijk](bc,a) = \sum(e) t_ki^ea (jb|ce) (4+)
                // W[ijk](bc,a) = \sum(e) J[j](bc,e) T[ki](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J2, T, 0,
                            (k * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](bc,a) -= \sum(m) t_jm^bc <ki|ma> (4-)
                // W[ijk](bc,a) -= \sum(m) T[j](m,bc) I[ki](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](ca,b) = \sum(e) t_ij^eb (kc|ae) (5+)
                // W[ijk](ca,b) = \sum(e) J[k](ca,e) T[ij](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J3, T, 0,
                            (i * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ca,b) -= \sum(m) t_km^ca <ij|mb> (5-)
                // W[ijk](ca,b) -= \sum(m) T[k](m,ca) I[ij](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA + b, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](cb,a) = \sum(e) t_ji^ea (kc|be) (6+)
                // W[ijk](cb,a) = \sum(e) J[k](cb,e) T[ji](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J3, T, 0,
                            (j * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](cb,a) -= \sum(m) t_km^cb <ji|ma> (6-)
                // W[ijk](cb,a) -= \sum(m) T[k](m,cb) I[ji](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA + a, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // V[ijk](ab,c) = W[ijk](ab,c)
                V->copy(W);

// V[ijk](ab,c) += t_i^a (jb|kc) + t_j^b (ia|kc) + t_k^c (ia|jb)
// V[ijk](ab,c) += f_ia T(jk|bc) + f_jb T(ik|ac) + f_kc T(ij|ab)
// Vt[ijk](ab,c) = V[ijk](ab,c) / (1 + \delta(abc))
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    long int ia = ia_idxAA->get(i, a);
                    for (long int b = 0; b < navirA; ++b) {
                        long int jb = ia_idxAA->get(j, b);
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int kc = ia_idxAA->get(k, c);
                            double value = V->get(ab, c) + (t1A->get(i, a) * J->get(jb, kc)) +
                                           (t1A->get(j, b) * J->get(ia, kc)) + (t1A->get(k, c) * J->get(ia, jb));

                            // E[4]_DT term
                            value += (FockA->get(i+nfrzc, a+noccA) * T->get(jk, bc)) +
                                           (FockA->get(j+nfrzc, b+noccA) * T->get(ik, ac)) + (FockA->get(k+nfrzc, c+noccA) * T->get(ij, ab));

                            double denom = 1 + ((a == b) + (b == c) + (a == c));
                            V->set(ab, c, value / denom);
                        }
                    }
                }

                // Denom
                double Dijk = Dij + FockA->get(k + nfrzc, k + nfrzc);
                double factor = 2 - ((i == j) + (j == k) + (i == k));

                // Compute energy
                double Xvalue, Yvalue, Zvalue;
#pragma omp parallel for private(Xvalue, Yvalue, Zvalue) reduction(+ : sum)
                for (long int a = 0; a < navirA; ++a) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b <= a; ++b) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c <= b; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);

                            // X_ijk^abc
                            Xvalue = (W->get(ab, c) * V->get(ab, c)) + (W->get(ac, b) * V->get(ac, b)) +
                                     (W->get(ba, c) * V->get(ba, c)) + (W->get(bc, a) * V->get(bc, a)) +
                                     (W->get(ca, b) * V->get(ca, b)) + (W->get(cb, a) * V->get(cb, a));

                            // Y_ijk^abc
                            Yvalue = V->get(ab, c) + V->get(bc, a) + V->get(ca, b);

                            // Z_ijk^abc
                            Zvalue = V->get(ac, b) + V->get(ba, c) + V->get(cb, a);

                            // contributions to energy
                            double value = (Yvalue - (2.0 * Zvalue)) * (W->get(ab, c) + W->get(bc, a) + W->get(ca, b));
                            value += (Zvalue - (2.0 * Yvalue)) * (W->get(ac, b) + W->get(ba, c) + W->get(cb, a));
                            value += 3.0 * Xvalue;
                            double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);
                            sum += (value * factor) / Dijkabc;
                        }
                    }
                }

                // progress counter
                ind += 1;
                double percent = static_cast<double>(ind) / static_cast<double>(Nijk) * 100.0;
                if (percent >= next_print) {
                    stop = std::time(nullptr);
                    next_print += step_print;
                    outfile->Printf("              %5.1lf  %8d s\n", percent,
                                    static_cast<int>(stop) - static_cast<int>(start));
                }

            }  // k
        }      // j
    }          // i
    T.reset();
    J.reset();
    W.reset();
    V.reset();
    J1.reset();
    J2.reset();
    J3.reset();
    K.reset();
    Jt.reset();
    L.reset();
    I.reset();

    // set energy
    E_t = sum;
    Eccsd_t = Eccsd + E_t;

}  // end ccsd_canonic_triples

//======================================================================
//       (T): incore
//======================================================================
void DFOCC::ccsd_canonic_triples_hm() {
    // defs
    SharedTensor2d K, L, M, I, J, T, U, Tau, W, X, Y, Z;
    SharedTensor2d V, J1, J2, J3, Jt;
    SharedTensor1d Eijk;
    long int Nijk;

    // progress counter
    std::time_t stop, start = std::time(nullptr);
    long int ind = 0;
    double step_print = 10.0;
    double next_print = step_print;

    // Find number of unique ijk combinations (i>=j>=k)
    Nijk = naoccA * (naoccA + 1) * (naoccA + 2) / 6;
    outfile->Printf("\tNumber of ijk combinations: %i \n", Nijk);

    // Memory: OV^3 + 2*O^2V^2 + 2*V^3 + O^3V + V^2N

    // Read t2 amps
    t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->sort(1324, t2, 1.0, 0.0);
    t2.reset();

    // Form (ij|ka)
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    M->read(psio_, PSIF_DFOCC_INTS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, K, M, 1.0, 0.0);
    K.reset();
    I = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Form (ia|jb)
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, M, M, 1.0, 0.0);

    // B(iaQ)
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ);
    L = M->transpose();
    M.reset();

    // B(Q,ab)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
    K->read(psio_, PSIF_DFOCC_INTS, true, true);

    // Form (ia|bc)
    J1 = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|BC)", naoccA, navirA, navirA, navirA);
    J1->gemm(false, false, L, K, 1.0, 0.0);
    K.reset();
    L.reset();

    // malloc W[ijk](abc)
    W = std::make_shared<Tensor2d>("W[IJK] <AB|C>", navirA * navirA, navirA);
    V = std::make_shared<Tensor2d>("V[IJK] <BA|C>", navirA * navirA, navirA);

    // B(Q,a>=b)
    // K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    // K->read(psio_, PSIF_DFOCC_INTS);

    // main loop
    E_t = 0.0;
    double sum = 0.0;
    for (long int i = 0; i < naoccA; ++i) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);
        for (long int j = 0; j <= i; ++j) {
            long int ij = ij_idxAA->get(i, j);
            double Dij = Di + FockA->get(j + nfrzc, j + nfrzc);
            for (long int k = 0; k <= j; ++k) {
                long int ik = ij_idxAA->get(i, k);
                long int jk = ij_idxAA->get(j, k);
                // W[ijk](ab,c) = \sum(e) t_jk^ec (ia|be) (1+)
                // W[ijk](ab,c) = \sum(e) J[i](ab,e) T[jk](ec)
                W->contract(false, false, navirA * navirA, navirA, navirA, J1, T, i * navirA * navirA * navirA,
                            (j * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ab,c) -= \sum(m) t_im^ab <jk|mc> (1-)
                // W[ijk](ab,c) -= \sum(m) T[i](m,ab) I[jk](mc)
                W->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);

                // W[ijk](ac,b) = \sum(e) t_kj^eb (ia|ce) (2+)
                // W[ijk](ac,b) = \sum(e) J[i](ac,e) T[kj](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, i * navirA * navirA * navirA,
                            (k * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ac,b) -= \sum(m) t_im^ac <kj|mb> (2-)
                // W[ijk](ac,b) -= \sum(m) T[i](m,ac) I[kj](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA * navirA + b, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](ba,c) = \sum(e) t_ik^ec (jb|ae) (3+)
                // W[ijk](ba,c) = \sum(e) J[j](ba,e) T[ik](ec)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, j * navirA * navirA * navirA,
                            (i * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ba,c) -= \sum(m) t_jm^ba <ik|mc> (3-)
                // W[ijk](ba,c) -= \sum(m) T[j](m,ba) I[ik](mc)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a * navirA, 1, V,
                                a * navirA * navirA + b * navirA, 1, 1.0);
                    }
                }

                // W[ijk](bc,a) = \sum(e) t_ki^ea (jb|ce) (4+)
                // W[ijk](bc,a) = \sum(e) J[j](bc,e) T[ki](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, j * navirA * navirA * navirA,
                            (k * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](bc,a) -= \sum(m) t_jm^bc <ki|ma> (4-)
                // W[ijk](bc,a) -= \sum(m) T[j](m,bc) I[ki](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](ca,b) = \sum(e) t_ij^eb (kc|ae) (5+)
                // W[ijk](ca,b) = \sum(e) J[k](ca,e) T[ij](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, k * navirA * navirA * navirA,
                            (i * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ca,b) -= \sum(m) t_km^ca <ij|mb> (5-)
                // W[ijk](ca,b) -= \sum(m) T[k](m,ca) I[ij](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA + b, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](cb,a) = \sum(e) t_ji^ea (kc|be) (6+)
                // W[ijk](cb,a) = \sum(e) J[k](cb,e) T[ji](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, k * navirA * navirA * navirA,
                            (j * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](cb,a) -= \sum(m) t_km^cb <ji|ma> (6-)
                // W[ijk](cb,a) -= \sum(m) T[k](m,cb) I[ji](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA + a, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // V[ijk](ab,c) = W[ijk](ab,c)
                V->copy(W);

// V[ijk](ab,c) += t_i^a (jb|kc) + t_j^b (ia|kc) + t_k^c (ia|jb)
// V[ijk](ab,c) += f_ia T(jk|bc) + f_jb T(ik|ac) + f_kc T(ij|ab)
// Vt[ijk](ab,c) = V[ijk](ab,c) / (1 + \delta(abc))
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    long int ia = ia_idxAA->get(i, a);
                    for (long int b = 0; b < navirA; ++b) {
                        long int jb = ia_idxAA->get(j, b);
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            long int kc = ia_idxAA->get(k, c);
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            double value = V->get(ab, c) + (t1A->get(i, a) * J->get(jb, kc)) +
                                           (t1A->get(j, b) * J->get(ia, kc)) + (t1A->get(k, c) * J->get(ia, jb));

                            // E[4]_DT term
                            value += (FockA->get(i+nfrzc, a+noccA) * T->get(jk, bc)) +
                                           (FockA->get(j+nfrzc, b+noccA) * T->get(ik, ac)) + (FockA->get(k+nfrzc, c+noccA) * T->get(ij, ab));

                            double denom = 1 + ((a == b) + (b == c) + (a == c));
                            V->set(ab, c, value / denom);
                        }
                    }
                }

                // Denom
                double Dijk = Dij + FockA->get(k + nfrzc, k + nfrzc);
                double factor = 2 - ((i == j) + (j == k) + (i == k));

                // Compute energy
                double Xvalue, Yvalue, Zvalue;
#pragma omp parallel for private(Xvalue, Yvalue, Zvalue) reduction(+ : sum)
                for (long int a = 0; a < navirA; ++a) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b <= a; ++b) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c <= b; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);

                            // X_ijk^abc
                            Xvalue = (W->get(ab, c) * V->get(ab, c)) + (W->get(ac, b) * V->get(ac, b)) +
                                     (W->get(ba, c) * V->get(ba, c)) + (W->get(bc, a) * V->get(bc, a)) +
                                     (W->get(ca, b) * V->get(ca, b)) + (W->get(cb, a) * V->get(cb, a));

                            // Y_ijk^abc
                            Yvalue = V->get(ab, c) + V->get(bc, a) + V->get(ca, b);

                            // Z_ijk^abc
                            Zvalue = V->get(ac, b) + V->get(ba, c) + V->get(cb, a);

                            // contributions to energy
                            double value = (Yvalue - (2.0 * Zvalue)) * (W->get(ab, c) + W->get(bc, a) + W->get(ca, b));
                            value += (Zvalue - (2.0 * Yvalue)) * (W->get(ac, b) + W->get(ba, c) + W->get(cb, a));
                            value += 3.0 * Xvalue;
                            double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);
                            sum += (value * factor) / Dijkabc;
                        }
                    }
                }

                // progress counter
                ind += 1;
                double percent = static_cast<double>(ind) / static_cast<double>(Nijk) * 100.0;
                if (percent >= next_print) {
                    stop = std::time(nullptr);
                    next_print += step_print;
                    outfile->Printf("              %5.1lf  %8d s\n", percent,
                                    static_cast<int>(stop) - static_cast<int>(start));
                }

            }  // k
        }      // j
    }          // i
    J1.reset();
    T.reset();
    J.reset();
    W.reset();
    V.reset();
    I.reset();

    // set energy
    E_t = sum;
    Eccsd_t = Eccsd + E_t;

}  // end ccsd_canonic_triples_hm

//======================================================================
//       (T): disk, This version includes E[4]_DT term
//======================================================================
void DFOCC::ccsd_canonic_triples_disk() {
    // defs
    SharedTensor2d K, L, M, I, J, T, U, Tau, W, X, Y, Z;
    SharedTensor2d V, J1, J2, J3, Jt;
    SharedTensor1d Eijk;
    long int Nijk;

    // progress counter
    std::time_t stop, start = std::time(nullptr);
    long int ind = 0;
    double step_print = 10.0;
    double next_print = step_print;

    // Find number of unique ijk combinations (i>=j>=k)
    Nijk = naoccA * (naoccA + 1) * (naoccA + 2) / 6;
    outfile->Printf("\tNumber of ijk combinations: %i \n", Nijk);

    // Malloc Eijk
    // Eijk = std::make_shared<Tensor1d>("Eijk", Nijk);

    // Memory: 2*O^2V^2 + 5*V^3 + O^3V + V^2N + V^3/2

    // Read t2 amps
    t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->sort(1324, t2, 1.0, 0.0);
    t2.reset();

    // Form (ij|ka)
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    M->read(psio_, PSIF_DFOCC_INTS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, K, M, 1.0, 0.0);
    K.reset();
    I = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Form (ia|jb)
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, M, M, 1.0, 0.0);

    // B(iaQ)
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ);
    L = M->transpose();
    M.reset();

    // malloc W[ijk](abc)
    W = std::make_shared<Tensor2d>("W[IJK] <AB|C>", navirA * navirA, navirA);
    V = std::make_shared<Tensor2d>("V[IJK] <BA|C>", navirA * navirA, navirA);
    J1 = std::make_shared<Tensor2d>("J[I] (A|BC)", navirA * navirA, navirA);
    J2 = std::make_shared<Tensor2d>("J[I] (A|BC)", navirA * navirA, navirA);
    J3 = std::make_shared<Tensor2d>("J[I] (A|BC)", navirA * navirA, navirA);

    // B(Q,ab)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    K->read(psio_, PSIF_DFOCC_INTS);

    // Form (ia|bc)
    Jt = std::make_shared<Tensor2d>("J[I] <A|B>=C", navirA, ntri_abAA);
    /*
    //psio_address addr = PSIO_ZERO;
    for(long int i = 0 ; i < naoccA; ++i){
        // Compute J[i](a,bc) = (ia|bc) = \sum(Q) B[i](aQ) * B(Q,bc)
        Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, i*navirA*nQ, 0, 1.0, 0.0);
        J1->expand23(navirA, navirA, navirA, Jt);

        // write
        psio_address addr = psio_get_address(PSIO_ZERO,(size_t)(i*navirA*navirA*navirA)*sizeof(double));
        J1->write(psio_, PSIF_DFOCC_INTS, addr, &addr);
    }
    */
    Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, 0, 0, 1.0, 0.0);
    J1->expand23(navirA, navirA, navirA, Jt);
    J1->mywrite(psio_, PSIF_DFOCC_IABC, false);
    for (long int i = 1; i < naoccA; ++i) {
        // Compute J[i](a,bc) = (ia|bc) = \sum(Q) B[i](aQ) * B(Q,bc)
        Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, i * navirA * nQ, 0, 1.0, 0.0);
        J1->expand23(navirA, navirA, navirA, Jt);

        // write
        J1->mywrite(psio_, PSIF_DFOCC_IABC, true);
    }
    K.reset();
    Jt.reset();
    L.reset();

    // main loop
    E_t = 0.0;
    double sum = 0.0;
    for (long int i = 0; i < naoccA; ++i) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);

        // Read J[i](a,bc)
        // psio_address addr1 = psio_get_address(PSIO_ZERO,(size_t)(i*navirA*navirA*navirA)*sizeof(double));
        // J1->read(psio_, PSIF_DFOCC_INTS, addr1, &addr1);
        J1->myread(psio_, PSIF_DFOCC_IABC, (size_t)(i * navirA * navirA * navirA) * sizeof(double));

        for (long int j = 0; j <= i; ++j) {
             long int ij = ij_idxAA->get(i, j);
             double Dij = Di + FockA->get(j + nfrzc, j + nfrzc);

            // Read J[j](a,bc)
            // psio_address addr2 = psio_get_address(PSIO_ZERO,(size_t)(j*navirA*navirA*navirA)*sizeof(double));
            // J2->read(psio_, PSIF_DFOCC_INTS, addr2, &addr2);
            J2->myread(psio_, PSIF_DFOCC_IABC, (size_t)(j * navirA * navirA * navirA) * sizeof(double));

            for (long int k = 0; k <= j; ++k) {
                 long int ik = ij_idxAA->get(i, k);
                 long int jk = ij_idxAA->get(j, k);
                // Read J[k](a,bc)
                // psio_address addr3 = psio_get_address(PSIO_ZERO,(size_t)(k*navirA*navirA*navirA)*sizeof(double));
                // J3->read(psio_, PSIF_DFOCC_INTS, addr3, &addr3);
                J3->myread(psio_, PSIF_DFOCC_IABC, (size_t)(k * navirA * navirA * navirA) * sizeof(double));

                // W[ijk](ab,c) = \sum(e) t_jk^ec (ia|be) (1+)
                // W[ijk](ab,c) = \sum(e) J[i](ab,e) T[jk](ec)
                W->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (j * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ab,c) -= \sum(m) t_im^ab <jk|mc> (1-)
                // W[ijk](ab,c) -= \sum(m) T[i](m,ab) I[jk](mc)
                W->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);

                // W[ijk](ac,b) = \sum(e) t_kj^eb (ia|ce) (2+)
                // W[ijk](ac,b) = \sum(e) J[i](ac,e) T[kj](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (k * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ac,b) -= \sum(m) t_im^ac <kj|mb> (2-)
                // W[ijk](ac,b) -= \sum(m) T[i](m,ac) I[kj](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA * navirA + b, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](ba,c) = \sum(e) t_ik^ec (jb|ae) (3+)
                // W[ijk](ba,c) = \sum(e) J[j](ba,e) T[ik](ec)
                V->contract(false, false, navirA * navirA, navirA, navirA, J2, T, 0,
                            (i * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ba,c) -= \sum(m) t_jm^ba <ik|mc> (3-)
                // W[ijk](ba,c) -= \sum(m) T[j](m,ba) I[ik](mc)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a * navirA, 1, V,
                                a * navirA * navirA + b * navirA, 1, 1.0);
                    }
                }

                // W[ijk](bc,a) = \sum(e) t_ki^ea (jb|ce) (4+)
                // W[ijk](bc,a) = \sum(e) J[j](bc,e) T[ki](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J2, T, 0,
                            (k * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](bc,a) -= \sum(m) t_jm^bc <ki|ma> (4-)
                // W[ijk](bc,a) -= \sum(m) T[j](m,bc) I[ki](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](ca,b) = \sum(e) t_ij^eb (kc|ae) (5+)
                // W[ijk](ca,b) = \sum(e) J[k](ca,e) T[ij](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J3, T, 0,
                            (i * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ca,b) -= \sum(m) t_km^ca <ij|mb> (5-)
                // W[ijk](ca,b) -= \sum(m) T[k](m,ca) I[ij](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA + b, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](cb,a) = \sum(e) t_ji^ea (kc|be) (6+)
                // W[ijk](cb,a) = \sum(e) J[k](cb,e) T[ji](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J3, T, 0,
                            (j * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](cb,a) -= \sum(m) t_km^cb <ji|ma> (6-)
                // W[ijk](cb,a) -= \sum(m) T[k](m,cb) I[ji](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA + a, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // V[ijk](ab,c) = W[ijk](ab,c)
                V->copy(W);

// V[ijk](ab,c) += t_i^a (jb|kc) + t_j^b (ia|kc) + t_k^c (ia|jb)
// V[ijk](ab,c) += f_ia T(jk|bc) + f_jb T(ik|ac) + f_kc T(ij|ab)
// Vt[ijk](ab,c) = V[ijk](ab,c) / (1 + \delta(abc))
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    long int ia = ia_idxAA->get(i, a);
                    for (long int b = 0; b < navirA; ++b) {
                        long int jb = ia_idxAA->get(j, b);
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            long int kc = ia_idxAA->get(k, c);
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            double value = V->get(ab, c) + (t1A->get(i, a) * J->get(jb, kc)) +
                                           (t1A->get(j, b) * J->get(ia, kc)) + (t1A->get(k, c) * J->get(ia, jb));

                            // E[4]_DT term
                            value += (FockA->get(i+nfrzc, a+noccA) * T->get(jk, bc)) +
                                           (FockA->get(j+nfrzc, b+noccA) * T->get(ik, ac)) + (FockA->get(k+nfrzc, c+noccA) * T->get(ij, ab));

                            double denom = 1 + ((a == b) + (b == c) + (a == c));
                            V->set(ab, c, value / denom);
                        }
                    }
                }

                // Denom
                double Dijk = Dij + FockA->get(k + nfrzc, k + nfrzc);
                double factor = 2 - ((i == j) + (j == k) + (i == k));

                // Compute energy
                double Xvalue, Yvalue, Zvalue;
#pragma omp parallel for private(Xvalue, Yvalue, Zvalue) reduction(+ : sum)
                for (long int a = 0; a < navirA; ++a) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b <= a; ++b) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c <= b; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);

                            // X_ijk^abc
                            Xvalue = (W->get(ab, c) * V->get(ab, c)) + (W->get(ac, b) * V->get(ac, b)) +
                                     (W->get(ba, c) * V->get(ba, c)) + (W->get(bc, a) * V->get(bc, a)) +
                                     (W->get(ca, b) * V->get(ca, b)) + (W->get(cb, a) * V->get(cb, a));

                            // Y_ijk^abc
                            Yvalue = V->get(ab, c) + V->get(bc, a) + V->get(ca, b);

                            // Z_ijk^abc
                            Zvalue = V->get(ac, b) + V->get(ba, c) + V->get(cb, a);

                            // contributions to energy
                            double value = (Yvalue - (2.0 * Zvalue)) * (W->get(ab, c) + W->get(bc, a) + W->get(ca, b));
                            value += (Zvalue - (2.0 * Yvalue)) * (W->get(ac, b) + W->get(ba, c) + W->get(cb, a));
                            value += 3.0 * Xvalue;
                            double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);
                            sum += (value * factor) / Dijkabc;
                        }
                    }
                }

                // progress counter
                ind += 1;
                double percent = static_cast<double>(ind) / static_cast<double>(Nijk) * 100.0;
                if (percent >= next_print) {
                    stop = std::time(nullptr);
                    next_print += step_print;
                    outfile->Printf("              %5.1lf  %8d s\n", percent,
                                    static_cast<int>(stop) - static_cast<int>(start));
                }

            }  // k
        }      // j
    }          // i
    T.reset();
    J.reset();
    W.reset();
    V.reset();
    J1.reset();
    J2.reset();
    J3.reset();
    I.reset();

    // set energy
    E_t = sum;
    Eccsd_t = Eccsd + E_t;

    // Delete the (IA|BC) file
    remove_binary_file(PSIF_DFOCC_IABC);

}  // end ccsd_canonic_triples_disk

//======================================================================
//       (T): grad
//======================================================================
void DFOCC::ccsd_canonic_triples_grad() {
    // defs
    SharedTensor2d K, L, M, I, I2, I3, J, T, U, Tau, W, X, Y, Z;
    SharedTensor2d V, J1, J2, J3, Jt, tL1, tL2, P2, P3, L2, L2c;
    SharedTensor2d Mijam, Mijab, Miabd;
    SharedTensor1d P1;
    long int Nijk;

    // Find number of unique ijk combinations (i>=j>=k)
    Nijk = naoccA * (naoccA + 1) * (naoccA + 2) / 6;
    outfile->Printf("\tNumber of ijk combinations: %i \n", Nijk);

    // Memory: 2*O^2V^2 + 4*V^3 + O^3V + V^2N + V^3/2

    // Read t2 amps
    if (!t2_incore) {
        t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    }
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->sort(1324, t2, 1.0, 0.0);
    t2.reset();

    // Form (ij|ka)
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    M->read(psio_, PSIF_DFOCC_INTS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, K, M, 1.0, 0.0);
    K.reset();
    I = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Form (ia|jb)
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, M, M, 1.0, 0.0);
    I2 = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA);
    I2->sort(1324, J, 1.0, 0.0);
    J.reset();

    // B(iaQ)
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ);
    L = M->transpose();
    M.reset();

    // B(Q,ab)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    K->read(psio_, PSIF_DFOCC_INTS);

    // Form (ia|bc)
    J1 = std::make_shared<Tensor2d>("J[I] (A|BC)", navirA * navirA, navirA);
    Jt = std::make_shared<Tensor2d>("J[I] <A|B>=C", navirA, ntri_abAA);
    Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, 0, 0, 1.0, 0.0);
    J1->expand23(navirA, navirA, navirA, Jt);
    J1->mywrite(psio_, PSIF_DFOCC_IABC, false);
    for (long int i = 1; i < naoccA; ++i) {
        // Compute J[i](a,bc) = (ia|bc) = \sum(Q) B[i](aQ) * B(Q,bc)
        Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, i * navirA * nQ, 0, 1.0, 0.0);
        J1->expand23(navirA, navirA, navirA, Jt);

        // write
        J1->mywrite(psio_, PSIF_DFOCC_IABC, true);
    }
    K.reset();
    Jt.reset();
    L.reset();

    // Alloc (t)^L_i^a amps
    tL1 = std::make_shared<Tensor2d>("(T)L <I|A>", naoccA, navirA);
    P1 = std::make_shared<Tensor1d>("P1 <A>", navirA);
    P2 = std::make_shared<Tensor2d>("P2 <A|B>", navirA, navirA);
    P3 = std::make_shared<Tensor2d>("P3 <A|I>", navirA, naoccA);
    // G1c_ii = std::make_shared<Tensor1d>("(T) Correlation OPDM <I|I>", naoccA);
    // G1c_aa = std::make_shared<Tensor1d>("(T) Correlation OPDM <A|A>", navirA);
    G1c_ia = std::make_shared<Tensor2d>("(T) Correlation OPDM <I|A>", naoccA, navirA);
    G1c_ab = std::make_shared<Tensor2d>("(T) Correlation OPDM <A|B>", navirA, navirA);
    L2 = std::make_shared<Tensor2d>("(T)AL2 <IJ|AB>", naoccA, naoccA, navirA, navirA);

    // malloc W[ijk](abc)
    W = std::make_shared<Tensor2d>("W[IJK] <AB|C>", navirA * navirA, navirA);
    V = std::make_shared<Tensor2d>("V[IJK] <AB|C>", navirA * navirA, navirA);
    I3 = std::make_shared<Tensor2d>("I[I] (A|BC)", navirA * navirA, navirA);
    Z = std::make_shared<Tensor2d>("Z[IJK] <I|AB>", naoccA, navirA * navirA);

    // Malloc M long intermediates
    Mijam = std::make_shared<Tensor2d>("M <IJ|AM>", naoccA, naoccA, navirA, naoccA);
    Mijab = std::make_shared<Tensor2d>("M <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Miabd = std::make_shared<Tensor2d>("M[I] <AB|D>", navirA * navirA, navirA);

    // main loop
    E_t = 0.0;
    double sum = 0.0;
    for (long int i = 0; i < naoccA; ++i) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);
        Miabd->zero();
        for (long int j = 0; j < naoccA; ++j) {
            long int ij = ij_idxAA->get(i, j);
            double Dij = Di + FockA->get(j + nfrzc, j + nfrzc);
            for (long int k = 0; k < naoccA; ++k) {
                long int ik = ij_idxAA->get(i, k);
                long int jk = ij_idxAA->get(j, k);

                // Read J[i](a,bc)
                J1->myread(psio_, PSIF_DFOCC_IABC, (size_t)(i * navirA * navirA * navirA) * sizeof(double));

                // W[ijk](ab,c) = \sum(e) t_jk^ec (ia|be) (1+):123
                // W[ijk](ab,c) = \sum(e) J[i](ab,e) T[jk](ec)
                W->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (j * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ab,c) -= \sum(m) t_im^ab <jk|mc> (1-)
                // W[ijk](ab,c) -= \sum(m) T[i](m,ab) I[jk](mc)
                W->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);

                // W[ijk](ac,b) = \sum(e) t_kj^eb (ia|ce) (2+):132
                // W[ijk](ac,b) = \sum(e) J[i](ac,e) T[kj](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (k * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ac,b) -= \sum(m) t_im^ac <kj|mb> (2-)
                // W[ijk](ac,b) -= \sum(m) T[i](m,ac) I[kj](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA * navirA + b, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // Read J[j](a,bc)
                J1->myread(psio_, PSIF_DFOCC_IABC, (size_t)(j * navirA * navirA * navirA) * sizeof(double));

                // W[ijk](ba,c) = \sum(e) t_ik^ec (jb|ae) (3+):213
                // W[ijk](ba,c) = \sum(e) J[j](ba,e) T[ik](ec)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (i * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ba,c) -= \sum(m) t_jm^ba <ik|mc> (3-)
                // W[ijk](ba,c) -= \sum(m) T[j](m,ba) I[ik](mc)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a * navirA, 1, V,
                                a * navirA * navirA + b * navirA, 1, 1.0);
                    }
                }

                // W[ijk](bc,a) = \sum(e) t_ki^ea (jb|ce) (4+):231
                // W[ijk](bc,a) = \sum(e) J[j](bc,e) T[ki](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (k * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](bc,a) -= \sum(m) t_jm^bc <ki|ma> (4-)
                // W[ijk](bc,a) -= \sum(m) T[j](m,bc) I[ki](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // Read J[k](a,bc)
                J1->myread(psio_, PSIF_DFOCC_IABC, (size_t)(k * navirA * navirA * navirA) * sizeof(double));

                // W[ijk](ca,b) = \sum(e) t_ij^eb (kc|ae) (5+):312
                // W[ijk](ca,b) = \sum(e) J[k](ca,e) T[ij](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (i * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ca,b) -= \sum(m) t_km^ca <ij|mb> (5-)
                // W[ijk](ca,b) -= \sum(m) T[k](m,ca) I[ij](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA + b, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](cb,a) = \sum(e) t_ji^ea (kc|be) (6+):321
                // W[ijk](cb,a) = \sum(e) J[k](cb,e) T[ji](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (j * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](cb,a) -= \sum(m) t_km^cb <ji|ma> (6-)
                // W[ijk](cb,a) -= \sum(m) T[k](m,cb) I[ji](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA + a, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // V[ijk](ab,c) = W[ijk](ab,c)
                V->copy(W);

// V[ijk](ab,c) += t_i^a (jb|kc) + t_j^b (ia|kc) + t_k^c (ia|jb)
// V[ijk](ab,c) += t_i^a <jk|bc> + t_j^b <ik|ac> + t_k^c <ij|ab>
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            double value = V->get(ab, c) + (t1A->get(i, a) * I2->get(jk, bc)) +
                                           (t1A->get(j, b) * I2->get(ik, ac)) + (t1A->get(k, c) * I2->get(ij, ab));
                            V->set(ab, c, value);
                        }
                    }
                }

                // Denom
                double Dijk = Dij + FockA->get(k + nfrzc, k + nfrzc);

// W = W/D and V = V/D
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b < navirA; ++b) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);
                            W->set(ab, c, W->get(ab, c) / Dijkabc);
                            V->set(ab, c, V->get(ab, c) / Dijkabc);
                        }
                    }
                }

                // Sort J[k](abc) -> I[k](bca)
                I3->sort3b(231, navirA, navirA, navirA, J1, 1.0, 0.0);

                // Compute energy
                double value_ = 0.0;
                double value2_ = 0.0;
                double value3_ = 0.0;
//#pragma omp parallel for private(value_,value2_,value3_) reduction(+:sum)
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b < navirA; ++b) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; ++c) {
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            long int ac = ab_idxAA->get(a, c);

                            // contributions to energy
                            value_ = (V->get(ab, c) - V->get(cb, a)) *
                                     ((4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b));
                            double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);
                            // sum += value_ / Dijkabc;
                            sum += value_ * Dijkabc;

                            // contributions to (t)^L_i^a
                            // (t)^L_i^a = 1/2 \sum(jkbc) X[ijk](a,bc) I[jk](bc)
                            value2_ = (4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) - (3.0 * W->get(cb, a)) -
                                      (2.0 * W->get(ac, b)) - W->get(ba, c);
                            // overwrite J1 with X
                            // J1->set(ab, c, value2_/Dijkabc);
                            J1->set(ab, c, value2_);

                            // contributions to G_ii & G_aa
                            // value3_ = V->get(ab,c) * ( (4.0*W->get(ab,c)) + W->get(bc,a) + W->get(ca,b) -
                            // (3.0*W->get(cb,a)) - (2.0*W->get(ac,b)) - W->get(ba,c) ) ;
                            /*
                            value3_ = value2_ * V->get(ab,c);
                            G1c_ii->subtract(i, value3_/(Dijkabc*Dijkabc));
                            G1c_aa->add(a, value3_/(Dijkabc*Dijkabc));
                            */
                        }
                    }
                }

                // Compute (t)^L_i^a
                P2->get_row(I2, jk);
                P1->gemv(false, navirA, navirA * navirA, J1, P2, 0.5, 0.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    tL1->add(i, a, P1->get(a));
                }

// G_ia = \sum(jkbc) X[ijk](a,bc) T[jk](bc)
// overwrite J1 with X
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            double value = (4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) - W->get(cb, a) -
                                           (2.0 * W->get(ac, b)) - (3.0 * W->get(ba, c));
                            J1->set(ab, c, value);
                        }
                    }
                }
                P2->get_row(T, jk);
                P1->gemv(false, navirA, navirA * navirA, J1, P2, 0.5, 0.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    G1c_ia->add(i, a, P1->get(a));
                }

// G_ae = -1/2 \sum(ijkbc) X[ijk](a,bc) Y[ijk](e,bc)
// overwrite J1 with X
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            double value = (4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) -
                                           (3.0 * W->get(cb, a)) - (2.0 * W->get(ac, b)) - W->get(ba, c);
                            J1->set(ab, c, value);
                        }
                    }
                }
                // G[ijk](a,e) += -1/2 \sum(bc) X[ijk](a,bc) Y[ijk](e,bc)
                G1c_ab->contract(false, true, navirA, navirA, navirA * navirA, J1, V, -0.5, 1.0);

// Mijab = \sum(kc) X[ijk](ab,c) T[k](c)
// overwrite J1 with X
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            double value = (4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) -
                                           (2.0 * W->get(cb, a)) - (2.0 * W->get(ac, b)) - (2.0 * W->get(ba, c));
                            J1->set(ab, c, value);
                        }
                    }
                }
                // M[ijk](a,b) = \sum(c) X[ijk](ab,c) T[k](c)
                P2->contract(false, true, navirA * navirA, 1, navirA, J1, t1A, 0, k * navirA, 1.0, 0.0);
                // M(ij,ab) = \sum(k) M[ijk](a,b)
                Mijab->add2row(P2, ij);

                //================================================
                //============ Form M_ijk^abc ====================
                //================================================
                // M[ijk](ab,c) = ( W[ijk](ab,c) + V[ijk](ab,c) ) / D_ijk^abc
                // Overwrite W with M
                W->axpy(V, 1.0);

// Compute (t)^L_ij^ab
// (t)^L_ij^ad += \sum(kbc) X[ijk](a,bc) I[k](d,bc)
// overwrite V with X
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            long int cb = ab_idxAA->get(c, b);
                            long int ac = ab_idxAA->get(a, c);
                            double value = (2.0 * W->get(ab, c)) - W->get(cb, a) - W->get(ac, b);
                            V->set(ab, c, value);
                        }
                    }
                }
                // P[ijk](a,d) = \sum(bc) X[ijk](a,bc) I[k](d,bc)
                P2->contract(false, true, navirA, navirA, navirA * navirA, V, I3, 1.0, 0.0);
                // (t)^L(ij,ad) += \sum(k) P[ijk](a,d)
                L2->add2row(P2, ij);

                // (t)^L_il^ab -= \sum(jkc) X[ijk](ab,c) I[jk](l,c)
                // Z[ijk](l,ab) = -\sum(c) I[jk](l,c) X[ijk](ab,c)
                Z->contract(false, true, naoccA, navirA * navirA, navirA, I, V,
                            (j * naoccA * naoccA * navirA) + (k * naoccA * navirA), 0, -1.0, 0.0);
                // (t)^L(il,ab) += \sum(jk) Z[ijk](l,ab)
                L2->axpy((size_t)naoccA * navirA * navirA, 0, 1, Z, i * naoccA * navirA * navirA, 1, 1.0);

// Mijam = \sum(kbc) X[ijk](a,cb) T[k](m,cb)
// overwrite V with X
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            double value = (2.0 * W->get(ab, c)) + (2.0 * W->get(bc, a)) + (2.0 * W->get(ca, b)) -
                                           W->get(cb, a) - W->get(ac, b) - (4.0 * W->get(ba, c));
                            V->set(ac, b, value);
                        }
                    }
                }
                // M[ijk](a,m) = \sum(bc) X[ijk](a,cb) T[k](m,cb)
                P3->contract(false, true, navirA, naoccA, navirA * navirA, V, T, 0, k * naoccA * navirA * navirA, 1.0,
                             0.0);
                // M(ij,am) = \sum(k) M[ijk](a,m)
                Mijam->add2row(P3, ij);

// overwrite V with X
// Miabd += \sum(jkc) X[ijk](ab,c) T[jk](d,c)
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            double value = (4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) -
                                           (2.0 * W->get(cb, a)) - (2.0 * W->get(ac, b)) - (2.0 * W->get(ba, c));
                            V->set(ab, c, value);
                        }
                    }
                }
                Miabd->contract(false, true, navirA * navirA, navirA, navirA, V, T, 0,
                                (j * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 1.0);

            }  // k
        }      // j
        Miabd->mywrite(psio_, PSIF_DFOCC_MIABC, true);
    }  // i
    // T.reset();
    I2.reset();
    W.reset();
    V.reset();
    J1.reset();
    I3.reset();
    Miabd.reset();
    I.reset();
    Z.reset();
    P1.reset();
    P2.reset();
    P3.reset();

    // set energy
    E_t = sum / 3.0;
    Eccsd_t = Eccsd + E_t;

    // Write M to disk
    Mijab->write(psio_, PSIF_DFOCC_DENS);
    Mijab.reset();
    Mijam->write(psio_, PSIF_DFOCC_DENS);
    Mijam.reset();

    // Form L2
    L2c = std::make_shared<Tensor2d>("(T)AL2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    L2c->sort(1324, L2, 1.0, 0.0);
    L2.reset();
    tL2 = std::make_shared<Tensor2d>("(T)L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    tL2->symmetrize(L2c);
    tL2->scale(2.0);
    L2c.reset();
    tL2->write_symm(psio_, PSIF_DFOCC_AMPS);
    tL2.reset();

    // write
    tL1->write(psio_, PSIF_DFOCC_AMPS);
    tL1.reset();

    // Delete the (IA|BC) file
    remove_binary_file(PSIF_DFOCC_IABC);

    //==========================================================================
    // Start for WIJK & VIJK
    //==========================================================================
    // Form (ij|ka)
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    M->read(psio_, PSIF_DFOCC_INTS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, K, M, 1.0, 0.0);
    K.reset();
    I = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <AI|JK>", navirA, naoccA, naoccA, naoccA);
    I->sort(4132, J, 1.0, 0.0);
    J.reset();

    // Form (ia|jb)
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, M, M, 1.0, 0.0);
    I2 = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA);
    I2->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Read
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AI)", nQ, navirA, naoccA);
    K->swap_3index_col(M);
    M.reset();

    // B(iaQ)
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (AI|Q)", naoccA * navirA, nQ);
    L = K->transpose();
    K.reset();

    // B(Q,ab)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    K->read(psio_, PSIF_DFOCC_INTS);

    // Form (ia|bc)
    J1 = std::make_shared<Tensor2d>("J[A] (IB|C)", naoccA * navirA, navirA);
    Jt = std::make_shared<Tensor2d>("J[A] <I|B>=C", naoccA, ntri_abAA);
    I3 = std::make_shared<Tensor2d>("I[AB] <I|C>", naoccA, navirA);
    for (long int a = 0; a < navirA; ++a) {
        // Compute J[a](i,bc) = (ai|bc) = \sum(Q) B[a](iQ) * B(Q,bc)
        Jt->contract(false, false, naoccA, ntri_abAA, nQ, L, K, a * naoccA * nQ, 0, 1.0, 0.0);
        J1->expand23(naoccA, navirA, navirA, Jt);

        //#pragma omp parallel for
        for (long int b = 0; b < navirA; ++b) {
            for (long int i = 0; i < naoccA; ++i) {
                long int ib = ia_idxAA->get(i, b);
                for (long int c = 0; c < navirA; ++c) {
                    I3->set(i, c, J1->get(ib, c));
                }
            }
            I3->mywrite(psio_, PSIF_DFOCC_ABIC, true);
        }
    }
    K.reset();
    Jt.reset();
    J1.reset();
    L.reset();

    // Form T transpose
    U = std::make_shared<Tensor2d>("T2 <AB|IJ>", navirA, navirA, naoccA, naoccA);
    U->trans(T);
    T.reset();

    // malloc W[abc](ijk)
    W = std::make_shared<Tensor2d>("W[ABC] <I|JK>", naoccA, naoccA * naoccA);
    V = std::make_shared<Tensor2d>("V[ABC] <I|JK>", naoccA, naoccA * naoccA);
    X = std::make_shared<Tensor2d>("X[ABC] <I|JK>", naoccA, naoccA * naoccA);
    G1c_ij = std::make_shared<Tensor2d>("(T) Correlation OPDM <I|J>", naoccA, naoccA);

    // main loop
    // E_t = 0.0;
    // sum = 0.0;
    size_t syc = 0;
    for (long int a = 0; a < navirA; ++a) {
        double Da = -1.0 * FockA->get(a + noccA, a + noccA);
        for (long int b = 0; b < navirA; ++b) {
            double Dab = Da - FockA->get(b + noccA, b + noccA);
            long int ab = ab_idxAA->get(a, b);
            for (long int c = 0; c < navirA; ++c) {
                long int ac = ab_idxAA->get(a, c);
                long int bc = ab_idxAA->get(b, c);

                // Read I[ac]
                syc = (a * navirA * naoccA * navirA) + (c * naoccA * navirA);
                I3->myread(psio_, PSIF_DFOCC_ABIC, syc * sizeof(double));

                // W[abc](i,jk) = \sum(e) I[ac](i,e) T'[b](e,jk) (1+):123
                W->contract(false, false, naoccA, naoccA * naoccA, navirA, I3, U, 0, b * navirA * naoccA * naoccA, 1.0,
                            0.0);

                // W[abc](i,jk) -= \sum(m) T'[ac](i,m) I[b](m,jk) (1-)
                W->contract(false, false, naoccA, naoccA * naoccA, naoccA, U, I,
                            (a * navirA * naoccA * naoccA) + (c * naoccA * naoccA), b * naoccA * naoccA * naoccA, -1.0,
                            1.0);

                // Read I[ab]
                syc = (a * navirA * naoccA * navirA) + (b * naoccA * navirA);
                I3->myread(psio_, PSIF_DFOCC_ABIC, syc * sizeof(double));

                // W[abc](i,kj) = \sum(e) I[ab](i,e) T'[c](e,kj) (2+):132
                V->contract(false, false, naoccA, naoccA * naoccA, navirA, I3, U, 0, c * navirA * naoccA * naoccA, 1.0,
                            0.0);

                // W[abc](i,kj) -= \sum(m) T'[ab](i,m) I[c](m,kj) (2-)
                V->contract(false, false, naoccA, naoccA * naoccA, naoccA, U, I,
                            (a * navirA * naoccA * naoccA) + (b * naoccA * naoccA), c * naoccA * naoccA * naoccA, -1.0,
                            1.0);
#pragma omp parallel for
                for (long int i = 0; i < naoccA; ++i) {
                    for (long int j = 0; j < naoccA; ++j) {
                        W->axpy((size_t)naoccA, i * naoccA * naoccA + j, naoccA, V, i * naoccA * naoccA + j * naoccA, 1,
                                1.0);
                    }
                }

                // Read I[bc]
                syc = (b * navirA * naoccA * navirA) + (c * naoccA * navirA);
                I3->myread(psio_, PSIF_DFOCC_ABIC, syc * sizeof(double));

                // W[abc](j,ik) = \sum(e) I[bc](j,e) T'[a](e,ik) (3+):213
                V->contract(false, false, naoccA, naoccA * naoccA, navirA, I3, U, 0, a * navirA * naoccA * naoccA, 1.0,
                            0.0);

                // W[abc](j,ik) -= \sum(m) T'[bc](j,m) I[a](m,ik) (3-)
                V->contract(false, false, naoccA, naoccA * naoccA, naoccA, U, I,
                            (b * navirA * naoccA * naoccA) + (c * naoccA * naoccA), a * naoccA * naoccA * naoccA, -1.0,
                            1.0);
#pragma omp parallel for
                for (long int i = 0; i < naoccA; ++i) {
                    for (long int j = 0; j < naoccA; ++j) {
                        W->axpy((size_t)naoccA, j * naoccA * naoccA + i * naoccA, 1, V,
                                i * naoccA * naoccA + j * naoccA, 1, 1.0);
                    }
                }

                // Read I[ba]
                syc = (b * navirA * naoccA * navirA) + (a * naoccA * navirA);
                I3->myread(psio_, PSIF_DFOCC_ABIC, syc * sizeof(double));

                // W[abc](j,ki) = \sum(e) I[ba](j,e) T'[c](e,ki) (4+):231
                V->contract(false, false, naoccA, naoccA * naoccA, navirA, I3, U, 0, c * navirA * naoccA * naoccA, 1.0,
                            0.0);

                // W[abc](j,ki) -= \sum(m) T'[ba](j,m) I[c](m,ki) (4-)
                V->contract(false, false, naoccA, naoccA * naoccA, naoccA, U, I,
                            (b * navirA * naoccA * naoccA) + (a * naoccA * naoccA), c * naoccA * naoccA * naoccA, -1.0,
                            1.0);
#pragma omp parallel for
                for (long int i = 0; i < naoccA; ++i) {
                    for (long int j = 0; j < naoccA; ++j) {
                        W->axpy((size_t)naoccA, j * naoccA * naoccA + i, naoccA, V, i * naoccA * naoccA + j * naoccA, 1,
                                1.0);
                    }
                }

                // Read I[cb]
                syc = (c * navirA * naoccA * navirA) + (b * naoccA * navirA);
                I3->myread(psio_, PSIF_DFOCC_ABIC, syc * sizeof(double));

                // W[abc](k,ij) = \sum(e) I[cb](k,e) T'[a](e,ij) (5+):312
                V->contract(false, false, naoccA, naoccA * naoccA, navirA, I3, U, 0, a * navirA * naoccA * naoccA, 1.0,
                            0.0);

                // W[abc](k,ij) -= \sum(m) T'[cb](k,m) I[a](m,ij) (5-)
                V->contract(false, false, naoccA, naoccA * naoccA, naoccA, U, I,
                            (c * navirA * naoccA * naoccA) + (b * naoccA * naoccA), a * naoccA * naoccA * naoccA, -1.0,
                            1.0);
#pragma omp parallel for
                for (long int i = 0; i < naoccA; ++i) {
                    for (long int j = 0; j < naoccA; ++j) {
                        W->axpy((size_t)naoccA, i * naoccA + j, naoccA * naoccA, V, i * naoccA * naoccA + j * naoccA, 1,
                                1.0);
                    }
                }

                // Read I[ca]
                syc = (c * navirA * naoccA * navirA) + (a * naoccA * navirA);
                I3->myread(psio_, PSIF_DFOCC_ABIC, syc * sizeof(double));

                // W[abc](k,ji) = \sum(e) I[ca](k,e) T'[b](e,ji) (6+):321
                V->contract(false, false, naoccA, naoccA * naoccA, navirA, I3, U, 0, b * navirA * naoccA * naoccA, 1.0,
                            0.0);

                // W[abc](k,ji) -= \sum(m) T'[ca](k,m) I[b](m,ji) (6-)
                V->contract(false, false, naoccA, naoccA * naoccA, naoccA, U, I,
                            (c * navirA * naoccA * naoccA) + (a * naoccA * naoccA), b * naoccA * naoccA * naoccA, -1.0,
                            1.0);
#pragma omp parallel for
                for (long int i = 0; i < naoccA; ++i) {
                    for (long int j = 0; j < naoccA; ++j) {
                        W->axpy((size_t)naoccA, j * naoccA + i, naoccA * naoccA, V, i * naoccA * naoccA + j * naoccA, 1,
                                1.0);
                    }
                }

                // V[abc](i,jk) = W[abc](i,jk)
                V->copy(W);

// V[abc](i,jk) += t_i^a (jb|kc) + t_j^b (ia|kc) + t_k^c (ia|jb)
// V[abc](i,jk) += t_i^a <jk|bc> + t_j^b <ik|ac> + t_k^c <ij|ab>
#pragma omp parallel for
                for (long int i = 0; i < naoccA; ++i) {
                    for (long int j = 0; j < naoccA; ++j) {
                        long int ij = ij_idxAA->get(i, j);
                        for (long int k = 0; k < naoccA; ++k) {
                            long int jk = ij_idxAA->get(j, k);
                            long int ik = ij_idxAA->get(i, k);
                            double value = V->get(i, jk) + (t1A->get(i, a) * I2->get(jk, bc)) +
                                           (t1A->get(j, b) * I2->get(ik, ac)) + (t1A->get(k, c) * I2->get(ij, ab));
                            V->set(i, jk, value);
                        }
                    }
                }

                // Denom
                double Dabc = Dab - FockA->get(c + noccA, c + noccA);

// W = W/D and V = V/D
#pragma omp parallel for
                for (long int i = 0; i < naoccA; ++i) {
                    double Diabc = Dabc + FockA->get(i + nfrzc, i + nfrzc);
                    for (long int j = 0; j < naoccA; ++j) {
                        double Dijabc = Diabc + FockA->get(j + nfrzc, j + nfrzc);
                        for (long int k = 0; k < naoccA; ++k) {
                            long int jk = ij_idxAA->get(j, k);
                            double Dijkabc = Dijabc + FockA->get(k + nfrzc, k + nfrzc);
                            W->set(i, jk, W->get(i, jk) / Dijkabc);
                            V->set(i, jk, V->get(i, jk) / Dijkabc);
                        }
                    }
                }

// Compute energy
/*
double value_ = 0.0;
//#pragma omp parallel for private(value_) reduction(+:sum)
#pragma omp parallel for
for(long int i = 0 ; i < naoccA; ++i){
    double Diabc = Dabc + FockA->get(i + nfrzc, i + nfrzc);
    for(long int j = 0 ; j < naoccA; ++j){
        double Dijabc = Diabc + FockA->get(j + nfrzc, j + nfrzc);
        long int ij = ij_idxAA->get(i,j);
        long int ji = ij_idxAA->get(j,i);
        for(long int k = 0 ; k < naoccA; ++k){
            double Dijkabc = Dijabc + FockA->get(k + nfrzc, k + nfrzc);
            long int jk = ij_idxAA->get(j,k);
            long int ki = ij_idxAA->get(k,i);

            // contributions to energy
            value_ = ( V->get(i,jk) - V->get(k,ji) ) * ( (4.0*W->get(i,jk)) + W->get(k,ij) + W->get(j,ki) ) ;
            sum += value_ * Dijkabc;
        }
    }
}
*/

// G_im = 1/2 \sum(jkabc) X[abc](i,jk) V[abc](m,jk)
// overwrite J1 with X
#pragma omp parallel for
                for (long int i = 0; i < naoccA; ++i) {
                    for (long int j = 0; j < naoccA; ++j) {
                        long int ij = ij_idxAA->get(i, j);
                        long int ji = ij_idxAA->get(j, i);
                        for (long int k = 0; k < naoccA; ++k) {
                            long int ik = ij_idxAA->get(i, k);
                            long int jk = ij_idxAA->get(j, k);
                            long int ki = ij_idxAA->get(k, i);
                            long int kj = ij_idxAA->get(k, j);
                            double value = (4.0 * W->get(i, jk)) + W->get(j, ki) + W->get(k, ij) -
                                           (3.0 * W->get(k, ji)) - (2.0 * W->get(i, kj)) - W->get(j, ik);
                            X->set(i, jk, value);
                        }
                    }
                }
                // G[abc](i,m) += 1/2 \sum(jk) X[abc](i,jk) V[abc](m,jk)
                G1c_ij->contract(false, true, naoccA, naoccA, naoccA * naoccA, X, V, 0.5, 1.0);

            }  // c
        }      // b
    }          // a
    U.reset();
    I2.reset();
    I3.reset();
    W.reset();
    V.reset();
    X.reset();
    I.reset();

    // set energy
    // E_t = sum/3.0;
    // Eccsd_t = Eccsd + E_t;

    // Delete the <AB|IC> file
    remove_binary_file(PSIF_DFOCC_ABIC);

    // Read t2 amps
    if (t2_incore || cc_lambda_ == "TRUE") {
        t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    }

}  // end ccsd_canonic_triples_grad

//======================================================================
//      New (T) Grad
//======================================================================
void DFOCC::ccsd_canonic_triples_grad2() {
    // defs
    SharedTensor2d K, L, M, I, I2, I3, J, T, U, Tau, W, X, Y, Z;
    SharedTensor2d V, J1, J2, J3, Jt, tL1, tL2, P2, P3, L2, L2c;
    SharedTensor2d Mijam, Mijab, Miabd;
    SharedTensor1d P1;

    long int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    outfile->Printf("\tnthreads: %i \n", nthreads);

    // This implementation requires the full set of ijk pairs
    long int Nijk = naoccA * naoccA * naoccA;
    outfile->Printf("\tNumber of ijk combinations: %i \n", Nijk);

    // Memory: 2*O^2V^2 + 4*V^3 + O^3V + V^2N + V^3/2

    // Read t2 amps
    if (!t2_incore) {
        t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    }
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->sort(1324, t2, 1.0, 0.0);
    t2.reset();

    // Form (ij|ka)
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    M->read(psio_, PSIF_DFOCC_INTS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, K, M, 1.0, 0.0);
    K.reset();
    I = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Form (ia|jb)
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, M, M, 1.0, 0.0);
    I2 = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|AB>", naoccA, naoccA, navirA, navirA);
    I2->sort(1324, J, 1.0, 0.0);
    J.reset();

    // B(iaQ)
    L = std::make_shared<Tensor2d>("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ);
    L = M->transpose();
    M.reset();

    // B(Q,ab)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    K->read(psio_, PSIF_DFOCC_INTS);

    // Form (ia|bc)
    J1 = std::make_shared<Tensor2d>("J[I] (A|BC)", navirA * navirA, navirA);
    Jt = std::make_shared<Tensor2d>("J[I] <A|B>=C", navirA, ntri_abAA);
    Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, 0, 0, 1.0, 0.0);
    J1->expand23(navirA, navirA, navirA, Jt);
    J1->mywrite(psio_, PSIF_DFOCC_IABC, false);
    for (long int i = 1; i < naoccA; ++i) {
        // Compute J[i](a,bc) = (ia|bc) = \sum(Q) B[i](aQ) * B(Q,bc)
        Jt->contract(false, false, navirA, ntri_abAA, nQ, L, K, i * navirA * nQ, 0, 1.0, 0.0);
        J1->expand23(navirA, navirA, navirA, Jt);

        // write
        J1->mywrite(psio_, PSIF_DFOCC_IABC, true);
    }
    K.reset();
    Jt.reset();
    L.reset();

    // Alloc (t)^L_i^a amps
    tL1 = std::make_shared<Tensor2d>("(T)L <I|A>", naoccA, navirA);
    P1 = std::make_shared<Tensor1d>("P1 <A>", navirA);
    P2 = std::make_shared<Tensor2d>("P2 <A|B>", navirA, navirA);
    P3 = std::make_shared<Tensor2d>("P3 <A|I>", navirA, naoccA);
    G1c_ii = std::make_shared<Tensor1d>("(T) Correlation OPDM <I|I>", naoccA);
    G1c_aa = std::make_shared<Tensor1d>("(T) Correlation OPDM <A|A>", navirA);
    L2 = std::make_shared<Tensor2d>("(T)AL2 <IJ|AB>", naoccA, naoccA, navirA, navirA);

    // malloc W[ijk](abc)
    W = std::make_shared<Tensor2d>("W[IJK] <AB|C>", navirA * navirA, navirA);
    V = std::make_shared<Tensor2d>("V[IJK] <AB|C>", navirA * navirA, navirA);
    I3 = std::make_shared<Tensor2d>("I[I] (A|BC)", navirA * navirA, navirA);
    Z = std::make_shared<Tensor2d>("Z[IJK] <I|AB>", naoccA, navirA * navirA);

    // Malloc M long intermediates
    Mijam = std::make_shared<Tensor2d>("M <IJ|AM>", naoccA, naoccA, navirA, naoccA);
    Mijab = std::make_shared<Tensor2d>("M <IJ|AB>", naoccA, naoccA, navirA, navirA);
    Miabd = std::make_shared<Tensor2d>("M[I] <AB|D>", navirA * navirA, navirA);

    // progress counter
    std::time_t stop, start = std::time(nullptr);
    long int ind = 0;
    double step_print = 10.0;
    double next_print = step_print;

    // main loop
    E_t = 0.0;
    double sum = 0.0;
    for (long int i = 0; i < naoccA; ++i) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);
        Miabd->zero();
        for (long int j = 0; j < naoccA; ++j) {
            long int ij = ij_idxAA->get(i, j);
            double Dij = Di + FockA->get(j + nfrzc, j + nfrzc);
            for (long int k = 0; k < naoccA; ++k) {
                long int ik = ij_idxAA->get(i, k);
                long int jk = ij_idxAA->get(j, k);

                // Read J[i](a,bc)
                J1->myread(psio_, PSIF_DFOCC_IABC, (size_t)(i * navirA * navirA * navirA) * sizeof(double));

                // W[ijk](ab,c) = \sum(e) t_jk^ec (ia|be) (1+):123
                // W[ijk](ab,c) = \sum(e) J[i](ab,e) T[jk](ec)
                W->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (j * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ab,c) -= \sum(m) t_im^ab <jk|mc> (1-)
                // W[ijk](ab,c) -= \sum(m) T[i](m,ab) I[jk](mc)
                W->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);

                // W[ijk](ac,b) = \sum(e) t_kj^eb (ia|ce) (2+):132
                // W[ijk](ac,b) = \sum(e) J[i](ac,e) T[kj](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (k * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ac,b) -= \sum(m) t_im^ac <kj|mb> (2-)
                // W[ijk](ac,b) -= \sum(m) T[i](m,ac) I[kj](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA * navirA + b, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // Read J[j](a,bc)
                J1->myread(psio_, PSIF_DFOCC_IABC, (size_t)(j * navirA * navirA * navirA) * sizeof(double));

                // W[ijk](ba,c) = \sum(e) t_ik^ec (jb|ae) (3+):213
                // W[ijk](ba,c) = \sum(e) J[j](ba,e) T[ik](ec)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (i * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ba,c) -= \sum(m) t_jm^ba <ik|mc> (3-)
                // W[ijk](ba,c) -= \sum(m) T[j](m,ba) I[ik](mc)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a * navirA, 1, V,
                                a * navirA * navirA + b * navirA, 1, 1.0);
                    }
                }

                // W[ijk](bc,a) = \sum(e) t_ki^ea (jb|ce) (4+):231
                // W[ijk](bc,a) = \sum(e) J[j](bc,e) T[ki](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (k * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](bc,a) -= \sum(m) t_jm^bc <ki|ma> (4-)
                // W[ijk](bc,a) -= \sum(m) T[j](m,bc) I[ki](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // Read J[k](a,bc)
                J1->myread(psio_, PSIF_DFOCC_IABC, (size_t)(k * navirA * navirA * navirA) * sizeof(double));

                // W[ijk](ca,b) = \sum(e) t_ij^eb (kc|ae) (5+):312
                // W[ijk](ca,b) = \sum(e) J[k](ca,e) T[ij](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (i * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ca,b) -= \sum(m) t_km^ca <ij|mb> (5-)
                // W[ijk](ca,b) -= \sum(m) T[k](m,ca) I[ij](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA + b, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](cb,a) = \sum(e) t_ji^ea (kc|be) (6+):321
                // W[ijk](cb,a) = \sum(e) J[k](cb,e) T[ji](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (j * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](cb,a) -= \sum(m) t_km^cb <ji|ma> (6-)
                // W[ijk](cb,a) -= \sum(m) T[k](m,cb) I[ji](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA + a, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // V[ijk](ab,c) = W[ijk](ab,c)
                V->copy(W);

// V[ijk](ab,c) += t_i^a (jb|kc) + t_j^b (ia|kc) + t_k^c (ia|jb)
// V[ijk](ab,c) += t_i^a <jk|bc> + t_j^b <ik|ac> + t_k^c <ij|ab>
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            double value = V->get(ab, c) + (t1A->get(i, a) * I2->get(jk, bc)) +
                                           (t1A->get(j, b) * I2->get(ik, ac)) + (t1A->get(k, c) * I2->get(ij, ab));
                            V->set(ab, c, value);
                        }
                    }
                }

                // Denom
                double Dijk = Dij + FockA->get(k + nfrzc, k + nfrzc);

// W = W/D and V = V/D
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b < navirA; ++b) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);
                            W->set(ab, c, W->get(ab, c) / Dijkabc);
                            V->set(ab, c, V->get(ab, c) / Dijkabc);
                        }
                    }
                }

                // Sort J[k](abc) -> I[k](bca)
                I3->sort3b(231, navirA, navirA, navirA, J1, 1.0, 0.0);

                // The following is better parallelized
                if (nthreads > 1) {
// Compute Energy & (t)^L_i^a
#pragma omp parallel for reduction(+ : sum)
                    for (long int a = 0; a < navirA; ++a) {
                        double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                        for (long int b = 0; b < navirA; ++b) {
                            double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                            long int ab = ab_idxAA->get(a, b);
                            long int ba = ab_idxAA->get(b, a);
                            for (long int c = 0; c < navirA; ++c) {
                                long int bc = ab_idxAA->get(b, c);
                                long int ca = ab_idxAA->get(c, a);
                                long int cb = ab_idxAA->get(c, b);
                                long int ac = ab_idxAA->get(a, c);

                                // contributions to energy
                                double value_ = (V->get(ab, c) - V->get(cb, a)) *
                                                ((4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b));
                                double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);
                                sum += value_ * Dijkabc;

                                // (t)^L_i^a = 1/2 \sum(jkbc) X[ijk](a,bc) I[jk](bc)
                                double value2_ = (4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) -
                                                 (3.0 * W->get(cb, a)) - (2.0 * W->get(ac, b)) - W->get(ba, c);
                                // overwrite J1 with X
                                J1->set(ab, c, value2_);
                            }
                        }
                    }

                    // Compute G_ii & G_aa
                    for (long int a = 0; a < navirA; ++a) {
                        double value3_ = 0.0;
#pragma omp parallel for reduction(+ : value3_)
                        for (long int b = 0; b < navirA; ++b) {
                            long int ab = ab_idxAA->get(a, b);
                            long int ba = ab_idxAA->get(b, a);
                            for (long int c = 0; c < navirA; ++c) {
                                long int bc = ab_idxAA->get(b, c);
                                long int ca = ab_idxAA->get(c, a);
                                long int cb = ab_idxAA->get(c, b);
                                long int ac = ab_idxAA->get(a, c);
                                value3_ +=
                                    V->get(ab, c) * ((4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) -
                                                     (3.0 * W->get(cb, a)) - (2.0 * W->get(ac, b)) - W->get(ba, c));
                            }
                        }
                        G1c_ii->subtract(i, value3_);
                        G1c_aa->add(a, value3_);
                    }
                }  // end if nthreads > 1

                // The following is better for serial jobs
                else {
                    // Compute energy
                    for (long int a = 0; a < navirA; ++a) {
                        double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                        for (long int b = 0; b < navirA; ++b) {
                            double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                            long int ab = ab_idxAA->get(a, b);
                            long int ba = ab_idxAA->get(b, a);
                            for (long int c = 0; c < navirA; ++c) {
                                long int bc = ab_idxAA->get(b, c);
                                long int ca = ab_idxAA->get(c, a);
                                long int cb = ab_idxAA->get(c, b);
                                long int ac = ab_idxAA->get(a, c);

                                // contributions to energy
                                double value_ = (V->get(ab, c) - V->get(cb, a)) *
                                                ((4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b));
                                double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);
                                // sum += value_ / Dijkabc;
                                sum += value_ * Dijkabc;

                                // contributions to (t)^L_i^a
                                // (t)^L_i^a = 1/2 \sum(jkbc) X[ijk](a,bc) I[jk](bc)
                                double value2_ = (4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) -
                                                 (3.0 * W->get(cb, a)) - (2.0 * W->get(ac, b)) - W->get(ba, c);

                                // overwrite J1 with X
                                // J1->set(ab, c, value2_/Dijkabc);
                                J1->set(ab, c, value2_);

                                // contributions to G_ii & G_aa
                                // value3_ = V->get(ab,c) * ( (4.0*W->get(ab,c)) + W->get(bc,a) + W->get(ca,b) -
                                // (3.0*W->get(cb,a)) - (2.0*W->get(ac,b)) - W->get(ba,c) ) ;
                                double value3_ = value2_ * V->get(ab, c);
                                G1c_ii->subtract(i, value3_);
                                G1c_aa->add(a, value3_);
                            }
                        }
                    }
                }  // end else

                // Compute (t)^L_i^a
                P2->get_row(I2, jk);
                P1->gemv(false, navirA, navirA * navirA, J1, P2, 0.5, 0.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    tL1->add(i, a, P1->get(a));
                }

// Mijab = \sum(kc) X[ijk](ab,c) T[k](c)
// overwrite J1 with X
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            double value = (4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) -
                                           (2.0 * W->get(cb, a)) - (2.0 * W->get(ac, b)) - (2.0 * W->get(ba, c));
                            J1->set(ab, c, value);
                        }
                    }
                }
                // M[ijk](a,b) = \sum(c) X[ijk](ab,c) T[k](c)
                P2->contract(false, true, navirA * navirA, 1, navirA, J1, t1A, 0, k * navirA, 1.0, 0.0);
                // M(ij,ab) = \sum(k) M[ijk](a,b)
                Mijab->add2row(P2, ij);

                //================================================
                //============ Form M_ijk^abc ====================
                //================================================
                // M[ijk](ab,c) = ( W[ijk](ab,c) + V[ijk](ab,c) ) / D_ijk^abc
                // Overwrite W with M
                W->axpy(V, 1.0);

// Compute (t)^L_ij^ab
// (t)^L_ij^ad += \sum(kbc) X[ijk](a,bc) I[k](d,bc)
// overwrite V with X
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            long int cb = ab_idxAA->get(c, b);
                            long int ac = ab_idxAA->get(a, c);
                            double value = (2.0 * W->get(ab, c)) - W->get(cb, a) - W->get(ac, b);
                            V->set(ab, c, value);
                        }
                    }
                }
                // P[ijk](a,d) = \sum(bc) X[ijk](a,bc) I[k](d,bc)
                P2->contract(false, true, navirA, navirA, navirA * navirA, V, I3, 1.0, 0.0);
                // (t)^L(ij,ad) += \sum(k) P[ijk](a,d)
                L2->add2row(P2, ij);

                // (t)^L_il^ab -= \sum(jkc) X[ijk](ab,c) I[jk](l,c)
                // Z[ijk](l,ab) = -\sum(c) I[jk](l,c) X[ijk](ab,c)
                Z->contract(false, true, naoccA, navirA * navirA, navirA, I, V,
                            (j * naoccA * naoccA * navirA) + (k * naoccA * navirA), 0, -1.0, 0.0);
                // (t)^L(il,ab) += \sum(jk) Z[ijk](l,ab)
                L2->axpy((size_t)naoccA * navirA * navirA, 0, 1, Z, i * naoccA * navirA * navirA, 1, 1.0);

// Mijam = \sum(kbc) X[ijk](a,cb) T[k](m,cb)
// overwrite V with X
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            double value = (2.0 * W->get(ab, c)) + (2.0 * W->get(bc, a)) + (2.0 * W->get(ca, b)) -
                                           W->get(cb, a) - W->get(ac, b) - (4.0 * W->get(ba, c));
                            V->set(ac, b, value);
                        }
                    }
                }
                // M[ijk](a,m) = \sum(bc) X[ijk](a,cb) T[k](m,cb)
                P3->contract(false, true, navirA, naoccA, navirA * navirA, V, T, 0, k * naoccA * navirA * navirA, 1.0,
                             0.0);
                // M(ij,am) = \sum(k) M[ijk](a,m)
                Mijam->add2row(P3, ij);

// overwrite V with X
// Miabd += \sum(jkc) X[ijk](ab,c) T[jk](d,c)
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c < navirA; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);
                            double value = (4.0 * W->get(ab, c)) + W->get(bc, a) + W->get(ca, b) -
                                           (2.0 * W->get(cb, a)) - (2.0 * W->get(ac, b)) - (2.0 * W->get(ba, c));
                            V->set(ab, c, value);
                        }
                    }
                }
                Miabd->contract(false, true, navirA * navirA, navirA, navirA, V, T, 0,
                                (j * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 1.0);

                // progress counter
                ind += 1;
                double percent = static_cast<double>(ind) / static_cast<double>(Nijk) * 100.0;
                if (percent >= next_print) {
                    stop = std::time(nullptr);
                    next_print += step_print;
                    outfile->Printf("              %5.1lf  %8d s\n", percent,
                                    static_cast<int>(stop) - static_cast<int>(start));
                }

            }  // k
        }      // j
        Miabd->mywrite(psio_, PSIF_DFOCC_MIABC, true);
    }  // i
    I2.reset();
    W.reset();
    V.reset();
    J1.reset();
    I3.reset();
    Miabd.reset();
    I.reset();
    Z.reset();
    P1.reset();
    P2.reset();
    P3.reset();

    // set energy
    E_t = sum / 3.0;
    Eccsd_t = Eccsd + E_t;

    // Write M to disk
    Mijab->write(psio_, PSIF_DFOCC_DENS);
    Mijab.reset();
    Mijam->write(psio_, PSIF_DFOCC_DENS);
    Mijam.reset();

    // Form L2
    L2c = std::make_shared<Tensor2d>("(T)AL2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    L2c->sort(1324, L2, 1.0, 0.0);
    L2.reset();
    tL2 = std::make_shared<Tensor2d>("(T)L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    tL2->symmetrize(L2c);
    tL2->scale(2.0);
    L2c.reset();

    // debug
    //L2c = std::make_shared<Tensor2d>("(T)L2 <Ij|Ab>", naoccA, naoccA, navirA, navirA);
    //L2c->sort(1324, tL2, 1.0, 0.0);
    //L2c->print();


    //auto L2aa = std::make_shared<Tensor2d>("(T)L2 <IJ||AB>", naoccA, naoccA, navirA, navirA);
    //L2aa->sort(2134, L2c, -1.0, 0.0);
    //L2aa->add(L2c);
    //L2aa->print();
    //L2c.reset();
    //L2aa.reset();
    // end debug

    tL2->write_symm(psio_, PSIF_DFOCC_AMPS);
    tL2.reset();

    // write
    tL1->write(psio_, PSIF_DFOCC_AMPS);
    //tL1->print();
    tL1.reset();

    // Delete the (IA|BC) file
    remove_binary_file(PSIF_DFOCC_IABC);

    // Read t2 amps
    if (t2_incore || cc_lambda_ == "TRUE") {
        t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    }

}  // end ccsd_canonic_triples_grad2

//======================================================================
//       Asymetric Triples: (AT)
//======================================================================
void DFOCC::ccsdl_canonic_triples_disk() {
    // defs
    SharedTensor2d K, L, M, I, J, T, U, Tau, W, WL, X, Y, Z;
    SharedTensor2d V, J1, J2, J3, Jt;
    SharedTensor1d Eijk;
    long int Nijk;

    // progress counter
    std::time_t stop, start = std::time(nullptr);
    long int ind = 0;
    double step_print = 10.0;
    double next_print = step_print;

    // Find number of unique ijk combinations (i>=j>=k)
    Nijk = naoccA * (naoccA + 1) * (naoccA + 2) / 6;
    outfile->Printf("\tNumber of ijk combinations: %i \n", Nijk);

    // Malloc Eijk
    // Eijk = std::make_shared<Tensor1d>("Eijk", Nijk);

    // Memory: 3*O^2V^2 + 5*V^3 + O^3V + V^2N + V^3/2

    // Read t2 amps
    t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    t2->read_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T->sort(1324, t2, 1.0, 0.0);
    t2.reset();

    // Read l2 amps
    l2 = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    l2->read_symm(psio_, PSIF_DFOCC_AMPS);
    L = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    L->sort(1324, l2, 1.0, 0.0);
    l2.reset();

    // Form (ij|ka)
    M = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    M->read(psio_, PSIF_DFOCC_INTS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->read(psio_, PSIF_DFOCC_INTS);
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|KA)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, K, M, 1.0, 0.0);
    K.reset();
    I = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <IJ|KA>", naoccA, naoccA, naoccA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Form (ia|jb)
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IA|JB)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, M, M, 1.0, 0.0);

    // B(iaQ)
    U = std::make_shared<Tensor2d>("DF_BASIS_CC B (IA|Q)", naoccA * navirA, nQ);
    U = M->transpose();
    M.reset();

    // malloc W[ijk](abc)
    W = std::make_shared<Tensor2d>("W[IJK] <AB|C>", navirA * navirA, navirA);
    WL = std::make_shared<Tensor2d>("WL[IJK] <AB|C>", navirA * navirA, navirA);
    V = std::make_shared<Tensor2d>("V[IJK] <BA|C>", navirA * navirA, navirA);
    J1 = std::make_shared<Tensor2d>("J[I] (A|BC)", navirA * navirA, navirA);
    J2 = std::make_shared<Tensor2d>("J[I] (A|BC)", navirA * navirA, navirA);
    J3 = std::make_shared<Tensor2d>("J[I] (A|BC)", navirA * navirA, navirA);

    // B(Q,ab)
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA);
    K->read(psio_, PSIF_DFOCC_INTS);

    // Form (ia|bc)
    Jt = std::make_shared<Tensor2d>("J[I] <A|B>=C", navirA, ntri_abAA);
    /*
    //psio_address addr = PSIO_ZERO;
    for(long int i = 0 ; i < naoccA; ++i){
        // Compute J[i](a,bc) = (ia|bc) = \sum(Q) B[i](aQ) * B(Q,bc)
        Jt->contract(false, false, navirA, ntri_abAA, nQ, U, K, i*navirA*nQ, 0, 1.0, 0.0);
        J1->expand23(navirA, navirA, navirA, Jt);

        // write
        psio_address addr = psio_get_address(PSIO_ZERO,(size_t)(i*navirA*navirA*navirA)*sizeof(double));
        J1->write(psio_, PSIF_DFOCC_INTS, addr, &addr);
    }
    */
    Jt->contract(false, false, navirA, ntri_abAA, nQ, U, K, 0, 0, 1.0, 0.0);
    J1->expand23(navirA, navirA, navirA, Jt);
    J1->mywrite(psio_, PSIF_DFOCC_IABC, false);
    for (long int i = 1; i < naoccA; ++i) {
        // Compute J[i](a,bc) = (ia|bc) = \sum(Q) B[i](aQ) * B(Q,bc)
        Jt->contract(false, false, navirA, ntri_abAA, nQ, U, K, i * navirA * nQ, 0, 1.0, 0.0);
        J1->expand23(navirA, navirA, navirA, Jt);

        // write
        J1->mywrite(psio_, PSIF_DFOCC_IABC, true);
    }

    K.reset();
    Jt.reset();
    U.reset();

    // main loop
    E_at = 0.0;
    double sum = 0.0;
    for (long int i = 0; i < naoccA; ++i) {
        double Di = FockA->get(i + nfrzc, i + nfrzc);

        // Read J[i](a,bc)
        // psio_address addr1 = psio_get_address(PSIO_ZERO,(size_t)(i*navirA*navirA*navirA)*sizeof(double));
        // J1->read(psio_, PSIF_DFOCC_INTS, addr1, &addr1);
        J1->myread(psio_, PSIF_DFOCC_IABC, (size_t)(i * navirA * navirA * navirA) * sizeof(double));

        for (long int j = 0; j <= i; ++j) {
            long int ij = ij_idxAA->get(i, j);
            double Dij = Di + FockA->get(j + nfrzc, j + nfrzc);

            // Read J[j](a,bc)
            // psio_address addr2 = psio_get_address(PSIO_ZERO,(size_t)(j*navirA*navirA*navirA)*sizeof(double));
            // J2->read(psio_, PSIF_DFOCC_INTS, addr2, &addr2);
            J2->myread(psio_, PSIF_DFOCC_IABC, (size_t)(j * navirA * navirA * navirA) * sizeof(double));

            for (long int k = 0; k <= j; ++k) {
                 long int ik = ij_idxAA->get(i, k);
                 long int jk = ij_idxAA->get(j, k);
                // Read J[k](a,bc)
                // psio_address addr3 = psio_get_address(PSIO_ZERO,(size_t)(k*navirA*navirA*navirA)*sizeof(double));
                // J3->read(psio_, PSIF_DFOCC_INTS, addr3, &addr3);
                J3->myread(psio_, PSIF_DFOCC_IABC, (size_t)(k * navirA * navirA * navirA) * sizeof(double));

                // W[ijk](ab,c) = \sum(e) t_jk^ec (ia|be) (1+)
                // W[ijk](ab,c) = \sum(e) J[i](ab,e) T[jk](ec)
                W->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (j * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ab,c) -= \sum(m) t_im^ab <jk|mc> (1-)
                // W[ijk](ab,c) -= \sum(m) T[i](m,ab) I[jk](mc)
                W->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);

                // W[ijk](ac,b) = \sum(e) t_kj^eb (ia|ce) (2+)
                // W[ijk](ac,b) = \sum(e) J[i](ac,e) T[kj](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, T, 0,
                            (k * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ac,b) -= \sum(m) t_im^ac <kj|mb> (2-)
                // W[ijk](ac,b) -= \sum(m) T[i](m,ac) I[kj](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, i * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA * navirA + b, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](ba,c) = \sum(e) t_ik^ec (jb|ae) (3+)
                // W[ijk](ba,c) = \sum(e) J[j](ba,e) T[ik](ec)
                V->contract(false, false, navirA * navirA, navirA, navirA, J2, T, 0,
                            (i * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ba,c) -= \sum(m) t_jm^ba <ik|mc> (3-)
                // W[ijk](ba,c) -= \sum(m) T[j](m,ba) I[ik](mc)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a * navirA, 1, V,
                                a * navirA * navirA + b * navirA, 1, 1.0);
                    }
                }

                // W[ijk](bc,a) = \sum(e) t_ki^ea (jb|ce) (4+)
                // W[ijk](bc,a) = \sum(e) J[j](bc,e) T[ki](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J2, T, 0,
                            (k * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](bc,a) -= \sum(m) t_jm^bc <ki|ma> (4-)
                // W[ijk](bc,a) -= \sum(m) T[j](m,bc) I[ki](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, j * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA * navirA + a, navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](ca,b) = \sum(e) t_ij^eb (kc|ae) (5+)
                // W[ijk](ca,b) = \sum(e) J[k](ca,e) T[ij](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J3, T, 0,
                            (i * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ca,b) -= \sum(m) t_km^ca <ij|mb> (5-)
                // W[ijk](ca,b) -= \sum(m) T[k](m,ca) I[ij](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, a * navirA + b, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                // W[ijk](cb,a) = \sum(e) t_ji^ea (kc|be) (6+)
                // W[ijk](cb,a) = \sum(e) J[k](cb,e) T[ji](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J3, T, 0,
                            (j * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](cb,a) -= \sum(m) t_km^cb <ji|ma> (6-)
                // W[ijk](cb,a) -= \sum(m) T[k](m,cb) I[ji](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, T, I, k * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        W->axpy((size_t)navirA, b * navirA + a, navirA * navirA, V, a * navirA * navirA + b * navirA, 1,
                                1.0);
                    }
                }

                //=========================
                // Asymmetric WL[ijk][abc]
                //=========================

                // W[ijk](ab,c) = \sum(e) l_jk^ec (ia|be) (1+)
                // W[ijk](ab,c) = \sum(e) J[i](ab,e) L[jk](ec)
                WL->contract(false, false, navirA * navirA, navirA, navirA, J1, L, 0,
                             (j * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ab,c) -= \sum(m) l_im^ab <jk|mc> (1-)
                // W[ijk](ab,c) -= \sum(m) L[i](m,ab) I[jk](mc)
                WL->contract(true, false, navirA * navirA, navirA, naoccA, L, I, i * naoccA * navirA * navirA,
                             (j * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);

                // W[ijk](ac,b) = \sum(e) l_kj^eb (ia|ce) (2+)
                // W[ijk](ac,b) = \sum(e) J[i](ac,e) L[kj](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J1, L, 0,
                            (k * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ac,b) -= \sum(m) l_im^ac <kj|mb> (2-)
                // W[ijk](ac,b) -= \sum(m) L[i](m,ac) I[kj](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, L, I, i * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        WL->axpy((size_t)navirA, a * navirA * navirA + b, navirA, V, a * navirA * navirA + b * navirA,
                                 1, 1.0);
                    }
                }

                // W[ijk](ba,c) = \sum(e) l_ik^ec (jb|ae) (3+)
                // W[ijk](ba,c) = \sum(e) J[j](ba,e) L[ik](ec)
                V->contract(false, false, navirA * navirA, navirA, navirA, J2, L, 0,
                            (i * naoccA * navirA * navirA) + (k * navirA * navirA), 1.0, 0.0);

                // W[ijk](ba,c) -= \sum(m) l_jm^ba <ik|mc> (3-)
                // W[ijk](ba,c) -= \sum(m) L[j](m,ba) I[ik](mc)
                V->contract(true, false, navirA * navirA, navirA, naoccA, L, I, j * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (k * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        WL->axpy((size_t)navirA, b * navirA * navirA + a * navirA, 1, V,
                                 a * navirA * navirA + b * navirA, 1, 1.0);
                    }
                }

                // W[ijk](bc,a) = \sum(e) l_ki^ea (jb|ce) (4+)
                // W[ijk](bc,a) = \sum(e) J[j](bc,e) L[ki](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J2, L, 0,
                            (k * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](bc,a) -= \sum(m) l_jm^bc <ki|ma> (4-)
                // W[ijk](bc,a) -= \sum(m) L[j](m,bc) I[ki](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, L, I, j * naoccA * navirA * navirA,
                            (k * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        WL->axpy((size_t)navirA, b * navirA * navirA + a, navirA, V, a * navirA * navirA + b * navirA,
                                 1, 1.0);
                    }
                }

                // W[ijk](ca,b) = \sum(e) l_ij^eb (kc|ae) (5+)
                // W[ijk](ca,b) = \sum(e) J[k](ca,e) L[ij](eb)
                V->contract(false, false, navirA * navirA, navirA, navirA, J3, L, 0,
                            (i * naoccA * navirA * navirA) + (j * navirA * navirA), 1.0, 0.0);

                // W[ijk](ca,b) -= \sum(m) l_km^ca <ij|mb> (5-)
                // W[ijk](ca,b) -= \sum(m) L[k](m,ca) I[ij](mb)
                V->contract(true, false, navirA * navirA, navirA, naoccA, L, I, k * naoccA * navirA * navirA,
                            (i * naoccA * naoccA * navirA) + (j * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        WL->axpy((size_t)navirA, a * navirA + b, navirA * navirA, V, a * navirA * navirA + b * navirA,
                                 1, 1.0);
                    }
                }

                // W[ijk](cb,a) = \sum(e) l_ji^ea (kc|be) (6+)
                // W[ijk](cb,a) = \sum(e) J[k](cb,e) L[ji](ea)
                V->contract(false, false, navirA * navirA, navirA, navirA, J3, L, 0,
                            (j * naoccA * navirA * navirA) + (i * navirA * navirA), 1.0, 0.0);

                // W[ijk](cb,a) -= \sum(m) l_km^cb <ji|ma> (6-)
                // W[ijk](cb,a) -= \sum(m) L[k](m,cb) I[ji](ma)
                V->contract(true, false, navirA * navirA, navirA, naoccA, L, I, k * naoccA * navirA * navirA,
                            (j * naoccA * naoccA * navirA) + (i * naoccA * navirA), -1.0, 1.0);
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    for (long int b = 0; b < navirA; ++b) {
                        WL->axpy((size_t)navirA, b * navirA + a, navirA * navirA, V, a * navirA * navirA + b * navirA,
                                 1, 1.0);
                    }
                }

                // V[ijk](ab,c) = WL[ijk](ab,c)
                V->copy(WL);

// V[ijk](ab,c) += l_i^a (jb|kc) + l_j^b (ia|kc) + l_k^c (ia|jb)
// Vt[ijk](ab,c) = V[ijk](ab,c) / (1 + \delta(abc))
#pragma omp parallel for
                for (long int a = 0; a < navirA; ++a) {
                    long int ia = ia_idxAA->get(i, a);
                    for (long int b = 0; b < navirA; ++b) {
                        long int jb = ia_idxAA->get(j, b);
                        long int ab = ab_idxAA->get(a, b);
                        for (long int c = 0; c < navirA; ++c) {
                            long int kc = ia_idxAA->get(k, c);
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            double value = V->get(ab, c) + (l1A->get(i, a) * J->get(jb, kc)) +
                                           (l1A->get(j, b) * J->get(ia, kc)) + (l1A->get(k, c) * J->get(ia, jb));

                            // E[4]_DT term
                            value += (FockA->get(i+nfrzc, a+noccA) * L->get(jk, bc)) +
                                           (FockA->get(j+nfrzc, b+noccA) * L->get(ik, ac)) + (FockA->get(k+nfrzc, c+noccA) * L->get(ij, ab));

                            double denom = 1 + ((a == b) + (b == c) + (a == c));
                            V->set(ab, c, value / denom);
                        }
                    }
                }

                // Denom
                double Dijk = Dij + FockA->get(k + nfrzc, k + nfrzc);
                double factor = 2 - ((i == j) + (j == k) + (i == k));

                // Compute energy
                double Xvalue, Yvalue, Zvalue;
#pragma omp parallel for private(Xvalue, Yvalue, Zvalue) reduction(+ : sum)
                for (long int a = 0; a < navirA; ++a) {
                    double Dijka = Dijk - FockA->get(a + noccA, a + noccA);
                    for (long int b = 0; b <= a; ++b) {
                        double Dijkab = Dijka - FockA->get(b + noccA, b + noccA);
                        long int ab = ab_idxAA->get(a, b);
                        long int ba = ab_idxAA->get(b, a);
                        for (long int c = 0; c <= b; ++c) {
                            long int ac = ab_idxAA->get(a, c);
                            long int bc = ab_idxAA->get(b, c);
                            long int ca = ab_idxAA->get(c, a);
                            long int cb = ab_idxAA->get(c, b);

                            // X_ijk^abc
                            Xvalue = (W->get(ab, c) * V->get(ab, c)) + (W->get(ac, b) * V->get(ac, b)) +
                                     (W->get(ba, c) * V->get(ba, c)) + (W->get(bc, a) * V->get(bc, a)) +
                                     (W->get(ca, b) * V->get(ca, b)) + (W->get(cb, a) * V->get(cb, a));

                            // Y_ijk^abc
                            Yvalue = V->get(ab, c) + V->get(bc, a) + V->get(ca, b);

                            // Z_ijk^abc
                            Zvalue = V->get(ac, b) + V->get(ba, c) + V->get(cb, a);

                            // contributions to energy
                            double value = (Yvalue - (2.0 * Zvalue)) * (W->get(ab, c) + W->get(bc, a) + W->get(ca, b));
                            value += (Zvalue - (2.0 * Yvalue)) * (W->get(ac, b) + W->get(ba, c) + W->get(cb, a));
                            value += 3.0 * Xvalue;
                            double Dijkabc = Dijkab - FockA->get(c + noccA, c + noccA);
                            sum += (value * factor) / Dijkabc;
                        }
                    }
                }
                // progress counter
                ind += 1;
                double percent = static_cast<double>(ind) / static_cast<double>(Nijk) * 100.0;
                if (percent >= next_print) {
                    stop = std::time(nullptr);
                    next_print += step_print;
                    outfile->Printf("              %5.1lf  %8d s\n", percent,
                                    static_cast<int>(stop) - static_cast<int>(start));
                }

            }  // k
        }      // j
    }          // i
    T.reset();
    L.reset();
    J.reset();
    W.reset();
    WL.reset();
    V.reset();
    J1.reset();
    J2.reset();
    J3.reset();
    I.reset();

    // set energy
    E_at = sum;
    Eccsd_at = Eccsd + E_at;

    // Delete the (IA|BC) file
    remove_binary_file(PSIF_DFOCC_IABC);

}  // end ccsdl_canonic_triples_disk

}  // namespace dfoccwave
}  // namespace psi
