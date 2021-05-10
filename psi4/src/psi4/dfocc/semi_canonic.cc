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

// Semicanonicalizing RHF Fock matrix by diagonalizing active-occupied (AOCC-AOCC) and active-virtual (AVIR-AVIR) blocks
#include "psi4/psi4-dec.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::semi_canonic() {
    // tell oter functions tat orbitals are already semi canonical.
    orbs_already_sc = 1;

    SharedTensor2d UooA = std::make_shared<Tensor2d>("UooA", naoccA, naoccA);
    SharedTensor2d UvvA = std::make_shared<Tensor2d>("UvvA", navirA, navirA);
    SharedTensor2d FockooA = std::make_shared<Tensor2d>("Fock <I|J>", naoccA, naoccA);
    SharedTensor2d FockvvA = std::make_shared<Tensor2d>("Fock <A|B>", navirA, navirA);

// Fockoo alpha spin case
#pragma omp parallel for
    for (int i = 0; i < naoccA; ++i) {
        for (int j = 0; j < naoccA; ++j) {
            FockooA->set(i, j, FockA->get(i + nfrzc, j + nfrzc));
        }
    }

// Fockvv alpha spin case
#pragma omp parallel for
    for (int a = 0; a < navirA; ++a) {
        for (int b = 0; b < navirA; ++b) {
            int aa = a + noccA;
            int bb = b + noccA;
            FockvvA->set(a, b, FockA->get(aa, bb));
        }
    }

    // Diagonalize Fock
    FockooA->diagonalize(UooA, eigooA, cutoff);
    FockvvA->diagonalize(UvvA, eigvvA, cutoff);

    // Print orbital energies
    if (occ_orb_energy == "TRUE" && mo_optimized == 1) {
        outfile->Printf("\n\n\tOCC Alpha Orbital Energies (a.u.) \n");
        outfile->Printf("\t  ---------------------------------- \n");

        // print occ orb energy
        outfile->Printf("\tAlpha occupied orbitals\n");
        for (int i = 0; i < naoccA; i++) {
            outfile->Printf("\t%2d %20.10f \n", i, eigooA->get(i));

        }  // end loop over naocc

        // print vir orb energy
        outfile->Printf("\n\tAlpha virtual orbitals\n");
        for (int i = 0; i < navirA; i++) {
            outfile->Printf("\t%2d %20.10f \n", i + noccA, eigvvA->get(i));

        }  // end loop over naocc

    }  // end main if

    // Build U
    UorbA->zero();

    // set to identity: it is necessary if we have frozen core or frozen virtual orbitals.
    UorbA->identity();

// Uoo contribution alpha spin case
#pragma omp parallel for
    for (int i = 0; i < naoccA; ++i) {
        for (int j = 0; j < naoccA; ++j) {
            UorbA->set(i + nfrzc, j + nfrzc, UooA->get(i, j));
        }
    }

// Uvv contribution alpha spin case
#pragma omp parallel for
    for (int a = 0; a < navirA; ++a) {
        for (int b = 0; b < navirA; ++b) {
            int aa = a + noccA;
            int bb = b + noccA;
            UorbA->set(aa, bb, UvvA->get(a, b));
        }
    }

    // Get new MOs
    SharedTensor2d Ca_new = std::make_shared<Tensor2d>("New alpha MO coefficients", nso_, nmo_);
    Ca_new->gemm(false, false, CmoA, UorbA, 1.0, 0.0);
    CmoA->copy(Ca_new);
    Ca_new.reset();

    // New fock
    auto FtempA = std::make_shared<Tensor2d>("temp", nmo_, nmo_);
    FtempA->copy(FockA);
    FockA->transform(FtempA, UorbA);
    //FtempA->print();
    //FockA->print();

    if (print_ > 2) {
        UorbA->print();
        CmoA->print();
    }

    UooA.reset();
    UvvA.reset();
    FockooA.reset();
    FockvvA.reset();

    //==========================================================================================
    //========================= UHF REFERENCE ==================================================
    //==========================================================================================
    if (reference_ == "UNRESTRICTED") {
        SharedTensor2d UooB = std::make_shared<Tensor2d>("UooB", naoccB, naoccB);
        SharedTensor2d UvvB = std::make_shared<Tensor2d>("UvvB", navirB, navirB);
        SharedTensor2d FockooB = std::make_shared<Tensor2d>("Fock <i|j>", naoccB, naoccB);
        SharedTensor2d FockvvB = std::make_shared<Tensor2d>("Fock <a|b>", navirB, navirB);

// Fockoo beta spin case
#pragma omp parallel for
        for (int i = 0; i < naoccB; ++i) {
            for (int j = 0; j < naoccB; ++j) {
                FockooB->set(i, j, FockB->get(i + nfrzc, j + nfrzc));
            }
        }

// Fockvv beta spin case
#pragma omp parallel for
        for (int a = 0; a < navirB; ++a) {
            for (int b = 0; b < navirB; ++b) {
                int aa = a + noccB;
                int bb = b + noccB;
                FockvvB->set(a, b, FockB->get(aa, bb));
            }
        }

        // Diagonalize Fock
        FockooB->diagonalize(UooB, eigooB, cutoff);
        FockvvB->diagonalize(UvvB, eigvvB, cutoff);

        // Print orbital energies
        if (occ_orb_energy == "TRUE" && mo_optimized == 1) {
            outfile->Printf("\n\n\tOCC Beta Orbital Energies (a.u.) \n");
            outfile->Printf("\t  ---------------------------------- \n");

            // print occ orb energy
            outfile->Printf("\tBeta occupied orbitals\n");
            for (int i = 0; i < naoccB; i++) {
                outfile->Printf("\t%2d %20.10f \n", i, eigooB->get(i));

            }  // end loop over naocc

            // print vir orb energy
            outfile->Printf("\n\tBeta virtual orbitals\n");
            for (int i = 0; i < navirB; i++) {
                outfile->Printf("\t%2d %20.10f \n", i + noccB, eigvvB->get(i));

            }  // end loop over naocc

        }  // end main if

        // Build U
        UorbB->zero();

        // set to identity: it is necessary if we have frozen core or frozen virtual orbitals.
        UorbB->identity();

// Uoo contribution beta spin case
#pragma omp parallel for
        for (int i = 0; i < naoccB; ++i) {
            for (int j = 0; j < naoccB; ++j) {
                UorbB->set(i + nfrzc, j + nfrzc, UooB->get(i, j));
            }
        }

// Uvv contribution beta spin case
#pragma omp parallel for
        for (int a = 0; a < navirB; ++a) {
            for (int b = 0; b < navirB; ++b) {
                int aa = a + noccB;
                int bb = b + noccB;
                UorbB->set(aa, bb, UvvB->get(a, b));
            }
        }

        // Get new MOs
        SharedTensor2d Cb_new = std::make_shared<Tensor2d>("New beta MO coefficients", nso_, nmo_);
        Cb_new->gemm(false, false, CmoB, UorbB, 1.0, 0.0);
        CmoB->copy(Cb_new);
        Cb_new.reset();

        // New fock
        auto FtempB = std::make_shared<Tensor2d>("temp", nmo_, nmo_);
        FtempB->copy(FockB);
        FockB->transform(FtempB, UorbB);

        if (print_ > 2) {
            UorbB->print();
            CmoB->print();
        }

        UooB.reset();
        UvvB.reset();
        FockooB.reset();
        FockvvB.reset();
    }  // end uhf

    // build mo coeff blocks
    mo_coeff_blocks();
}//semi-canonic

//======================================================================
//       T2-trans
//======================================================================
void DFOCC::t2_trans(std::string amp_type) {
    SharedTensor2d T, X, Y, Z;

    SharedTensor2d UooA = std::make_shared<Tensor2d>("UooA", naoccA, naoccA);
    SharedTensor2d UvvA = std::make_shared<Tensor2d>("UvvA", navirA, navirA);

// Uoo contribution alpha spin case
#pragma omp parallel for
    for (int i = 0; i < naoccA; ++i) {
        for (int j = 0; j < naoccA; ++j) {
            UooA->set(i, j, UorbA->get(i + nfrzc, j + nfrzc));
        }
    }

// Uvv contribution alpha spin case
#pragma omp parallel for
    for (int a = 0; a < navirA; ++a) {
        for (int b = 0; b < navirA; ++b) {
            int aa = a + noccA;
            int bb = b + noccA;
            UvvA->set(a, b, UorbA->get(aa, bb));
        }
    }


    //RHF
    if (reference_ == "RESTRICTED") {
        // Read t2 amps
        if (amp_type == "T2") {
            t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        }
        else if (amp_type == "L2") {
            t2 = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        }
        t2->read_symm(psio_, PSIF_DFOCC_AMPS);
        T = std::make_shared<Tensor2d>("Non-canonic T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T->sort(1324, t2, 1.0, 0.0);
        t2.reset();

        // X(kl,cb) = \sum(d) T(kl,cd) * U(d,b)
        X = std::make_shared<Tensor2d>("Nc-sc T2 <KL|CB>", naoccA, naoccA, navirA, navirA);
        X->contract(false, false, naoccA * naoccA * navirA, navirA, navirA, T, UvvA, 1.0, 0.0);
        T.reset();

        // Y(kl,ab) = \sum(c) X(kl,cb) * U(c,a)
        T = std::make_shared<Tensor2d>("Nc-sc T2 <KL|BC>", naoccA, naoccA, navirA, navirA);
        T->sort(1243, X, 1.0, 0.0);
        X.reset();
        X = std::make_shared<Tensor2d>("Nc-sc T2 <KL|BA>", naoccA, naoccA, navirA, navirA);
        X->contract(false, false, naoccA * naoccA * navirA, navirA, navirA, T, UvvA, 1.0, 0.0);
        T.reset();
        Y = std::make_shared<Tensor2d>("Nc-sc T2 <AB|KL>", navirA, navirA, naoccA, naoccA);
        Y->sort(4312, X, 1.0, 0.0);
        X.reset();

        // Z(ab,kj) = \sum(l) Y(ab,kl) * U(l,j)
        Z = std::make_shared<Tensor2d>("Nc-sc T2 <AB|KJ>", navirA, navirA, naoccA, naoccA);
        Z->contract(false, false, navirA * navirA * naoccA, naoccA, naoccA, Y, UooA, 1.0, 0.0);
        Y.reset();

        // T(ab,ij) = \sum(k) Z(ab,kj) * U(k,i)
        T = std::make_shared<Tensor2d>("Nc-sc T2 <AB|JK>", navirA, navirA, naoccA, naoccA);
        T->sort(1243, Z, 1.0, 0.0);
        Z.reset();
        X = std::make_shared<Tensor2d>("Nc-sc T2 <AB|JI>", navirA, navirA, naoccA, naoccA);
        X->contract(false, false, navirA * navirA * naoccA, naoccA, naoccA, T, UooA, 1.0, 0.0);
        T.reset();
        T = std::make_shared<Tensor2d>("Semi-canonic T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        T->sort(4312, X, 1.0, 0.0);
        X.reset();

        // Now sort it to chemist notation since rhf code expect that
        if (amp_type == "T2") {
            t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        }
        else if (amp_type == "L2") {
            t2 = std::make_shared<Tensor2d>("L2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        }
        t2->sort(1324, T, 1.0, 0.0);
        T.reset();
        t2->write_symm(psio_, PSIF_DFOCC_AMPS);
        t2.reset();

    }//rhf
    //UHF
    else if (reference_ == "UNRESTRICTED") {

        SharedTensor2d UooB = std::make_shared<Tensor2d>("UooB", naoccB, naoccB);
        SharedTensor2d UvvB = std::make_shared<Tensor2d>("UvvB", navirB, navirB);

        // Uoo contribution beta spin case
#pragma omp parallel for
        for (int i = 0; i < naoccB; ++i) {
            for (int j = 0; j < naoccB; ++j) {
                UooB->set(i, j, UorbB->get(i + nfrzc, j + nfrzc));
            }
        }

        // Uvv contribution beta spin case
#pragma omp parallel for
        for (int a = 0; a < navirB; ++a) {
            for (int b = 0; b < navirB; ++b) {
                int aa = a + noccB;
                int bb = b + noccB;
                UvvB->set(a, b, UorbB->get(aa, bb));
            }
        }

        // T2AA
        // Read t2 amps
        if (amp_type == "T2") {
            T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        }
        else if (amp_type == "L2") {
            T = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        }
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

        // X(kl,cb) = \sum(d) T(kl,cd) * U(d,b)
        X = std::make_shared<Tensor2d>("Nc-sc T2 <KL|CB>", naoccA, naoccA, navirA, navirA);
        X->contract(false, false, naoccA * naoccA * navirA, navirA, navirA, T, UvvA, 1.0, 0.0);
        T.reset();

        // Y(kl,ab) = \sum(c) X(kl,cb) * U(c,a)
        T = std::make_shared<Tensor2d>("Nc-sc T2 <KL|BC>", naoccA, naoccA, navirA, navirA);
        T->sort(1243, X, 1.0, 0.0);
        X.reset();
        X = std::make_shared<Tensor2d>("Nc-sc T2 <KL|BA>", naoccA, naoccA, navirA, navirA);
        X->contract(false, false, naoccA * naoccA * navirA, navirA, navirA, T, UvvA, 1.0, 0.0);
        T.reset();
        Y = std::make_shared<Tensor2d>("Nc-sc T2 <AB|KL>", navirA, navirA, naoccA, naoccA);
        Y->sort(4312, X, 1.0, 0.0);
        X.reset();

        // Z(ab,kj) = \sum(l) Y(ab,kl) * U(l,j)
        Z = std::make_shared<Tensor2d>("Nc-sc T2 <AB|KJ>", navirA, navirA, naoccA, naoccA);
        Z->contract(false, false, navirA * navirA * naoccA, naoccA, naoccA, Y, UooA, 1.0, 0.0);
        Y.reset();

        // T(ab,ij) = \sum(k) Z(ab,kj) * U(k,i)
        T = std::make_shared<Tensor2d>("Nc-sc T2 <AB|JK>", navirA, navirA, naoccA, naoccA);
        T->sort(1243, Z, 1.0, 0.0);
        Z.reset();
        X = std::make_shared<Tensor2d>("Nc-sc T2 <AB|JI>", navirA, navirA, naoccA, naoccA);
        X->contract(false, false, navirA * navirA * naoccA, naoccA, naoccA, T, UooA, 1.0, 0.0);
        T.reset();
        if (amp_type == "T2") {
            T = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        }
        else if (amp_type == "L2") {
            T = std::make_shared<Tensor2d>("L2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
        }
        T->sort(4312, X, 1.0, 0.0);
        X.reset();

        T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // T2BB
        // Read t2 amps
        if (amp_type == "T2") {
            T = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        }
        else if (amp_type == "L2") {
            T = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        }
        T->read_anti_symm(psio_, PSIF_DFOCC_AMPS);

        // X(kl,cb) = \sum(d) T(kl,cd) * U(d,b)
        X = std::make_shared<Tensor2d>("Nc-sc T2 <kl|cb>", naoccB, naoccB, navirB, navirB);
        X->contract(false, false, naoccB * naoccB * navirB, navirB, navirB, T, UvvB, 1.0, 0.0);
        T.reset();

        // Y(kl,ab) = \sum(c) X(kl,cb) * U(c,a)
        T = std::make_shared<Tensor2d>("Nc-sc T2 <kl|bc>", naoccB, naoccB, navirB, navirB);
        T->sort(1243, X, 1.0, 0.0);
        X.reset();
        X = std::make_shared<Tensor2d>("Nc-sc T2 <kl|ba>", naoccB, naoccB, navirB, navirB);
        X->contract(false, false, naoccB * naoccB * navirB, navirB, navirB, T, UvvB, 1.0, 0.0);
        T.reset();
        Y = std::make_shared<Tensor2d>("Nc-sc T2 <ab|kl>", navirB, navirB, naoccB, naoccB);
        Y->sort(4312, X, 1.0, 0.0);
        X.reset();

        // Z(ab,kj) = \sum(l) Y(ab,kl) * U(l,j)
        Z = std::make_shared<Tensor2d>("Nc-sc T2 <ab|kj>", navirB, navirB, naoccB, naoccB);
        Z->contract(false, false, navirB * navirB * naoccB, naoccB, naoccB, Y, UooB, 1.0, 0.0);
        Y.reset();

        // T(ab,ij) = \sum(k) Z(ab,kj) * U(k,i)
        T = std::make_shared<Tensor2d>("Nc-sc T2 <ab|jk>", navirB, navirB, naoccB, naoccB);
        T->sort(1243, Z, 1.0, 0.0);
        Z.reset();
        X = std::make_shared<Tensor2d>("Nc-sc T2 <ab|ji>", navirB, navirB, naoccB, naoccB);
        X->contract(false, false, navirB * navirB * naoccB, naoccB, naoccB, T, UooB, 1.0, 0.0);
        T.reset();
        if (amp_type == "T2") {
            T = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        }
        else if (amp_type == "L2") {
            T = std::make_shared<Tensor2d>("L2 <ij|ab>", naoccB, naoccB, navirB, navirB);
        }
        T->sort(4312, X, 1.0, 0.0);
        X.reset();

        T->write_anti_symm(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        // T2AB
        // Read t2 amps
        if (amp_type == "T2") {
            T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        }
        else if (amp_type == "L2") {
            T = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        }
        T->read(psio_, PSIF_DFOCC_AMPS);

        // X(kl,cb) = \sum(d) T(kl,cd) * U(d,b)
        X = std::make_shared<Tensor2d>("Nc-sc T2 <Kl|Cb>", naoccA, naoccB, navirA, navirB);
        X->contract(false, false, naoccA * naoccB * navirA, navirB, navirB, T, UvvB, 1.0, 0.0);
        T.reset();

        // Y(kl,ab) = \sum(c) X(kl,cb) * U(c,a)
        T = std::make_shared<Tensor2d>("Nc-sc T2 <Kl|bC>", naoccA, naoccB, navirB, navirA);
        T->sort(1243, X, 1.0, 0.0);
        X.reset();
        X = std::make_shared<Tensor2d>("Nc-sc T2 <Kl|bA>", naoccA, naoccB, navirB, navirA);
        X->contract(false, false, naoccA * naoccB * navirB, navirA, navirA, T, UvvA, 1.0, 0.0);
        T.reset();
        Y = std::make_shared<Tensor2d>("Nc-sc T2 <Ab|Kl>", navirA, navirB, naoccA, naoccB);
        Y->sort(4312, X, 1.0, 0.0);
        X.reset();

        // Z(ab,kj) = \sum(l) Y(ab,kl) * U(l,j)
        Z = std::make_shared<Tensor2d>("Nc-sc T2 <Ab|Kj>", navirA, navirB, naoccA, naoccB);
        Z->contract(false, false, navirA * navirB * naoccA, naoccB, naoccB, Y, UooB, 1.0, 0.0);
        Y.reset();

        // T(ab,ij) = \sum(k) Z(ab,kj) * U(k,i)
        T = std::make_shared<Tensor2d>("Nc-sc T2 <Ab|jK>", navirA, navirB, naoccB, naoccA);
        T->sort(1243, Z, 1.0, 0.0);
        Z.reset();
        X = std::make_shared<Tensor2d>("Nc-sc T2 <Ab|jI>", navirA, navirB, naoccB, naoccA);
        X->contract(false, false, navirA * navirB * naoccB, naoccA, naoccA, T, UooA, 1.0, 0.0);
        T.reset();
        if (amp_type == "T2") {
            T = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        }
        else if (amp_type == "L2") {
            T = std::make_shared<Tensor2d>("L2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
        }
        T->sort(4312, X, 1.0, 0.0);
        X.reset();

        T->write(psio_, PSIF_DFOCC_AMPS);
        T.reset();

        UooB.reset();
        UvvB.reset();

    }//uhf

    UooA.reset();
    UvvA.reset();

}//t2_trans() {

}  // namespace dfoccwave
}  // namespace psi
