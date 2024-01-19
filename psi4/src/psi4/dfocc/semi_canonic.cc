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
}
}  // namespace dfoccwave
}  // namespace psi
