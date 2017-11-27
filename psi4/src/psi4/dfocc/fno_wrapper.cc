/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"
#include "fno.h"
#include "psi4/libmints/matrix.h"
#include "psi4/physconst.h"
#include "psi4/liboptions/liboptions.h"

namespace psi {
namespace dfoccwave {

void DFOCC::fno_wrapper() {
    outfile->Printf("\n\tComputing FNOs...\n");

    SharedTensor2d T;
    SharedTensor1d diag;

    // Trans b_ia^Q
    // Read SO integrals
    bQso = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_);
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);

    // Form B(Q,ia)
    if (do_cd == "FALSE")
        b_ia();
    else
        b_ia_cd();
    bQso.reset();

    // Read Ints
    bQiaA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
    bQiaA->read(psio_, PSIF_DFOCC_INTS);

    // Remove old integrals file
    remove_binary_file(PSIF_DFOCC_INTS);
    // outfile->Printf("\tI am here\n");

    timer_on("fno");
    if (reference_ == "RESTRICTED") {
        // malloc
        // T = std::make_shared<Tensor2d>("FNO T matrix <V|V>", nvirA, nvirA);
        TfnoA = std::make_shared<Tensor2d>("FNO T matrix", nmo_, nmo_);
        VfnoA = std::make_shared<Tensor2d>("FNO V matrix", nso_, nmo_);
        diag = std::make_shared<Tensor1d>("Occupation numbers", nvirA);

        // Call FnoCC
        SharedFnoCC fnoA;
        int navir_fno = 0;
        if (options_["ACTIVE_NAT_ORBS"].has_changed()) {
            for (int h = 0; h < nirrep_; h++) {
                navir_fno += (int)options_["ACTIVE_NAT_ORBS"][h].to_double();
            }
            fnoA = std::make_shared<FnoCC>("RHF FNO", nQ, nso_, nfrzc, naoccA, nvirA, CmoA, FockA, bQiaA, tol_fno,
                                         fno_percentage, navir_fno, false, true);
        }

        else if (options_["OCC_PERCENTAGE"].has_changed()) {
            fnoA = std::make_shared<FnoCC>("RHF FNO", nQ, nso_, nfrzc, naoccA, nvirA, CmoA, FockA, bQiaA, tol_fno,
                                         fno_percentage, navir_fno, true, false);
        }

        else {
            fnoA = std::make_shared<FnoCC>("RHF FNO", nQ, nso_, nfrzc, naoccA, nvirA, CmoA, FockA, bQiaA, tol_fno,
                                         fno_percentage, navir_fno, false, false);
        }

        // Get MP2 energy
        Ecorr = fnoA->mp2_corr();
        Emp2L = Eref + Ecorr;
        outfile->Printf("\tComputing DF-MP2 energy in the full virtual space... \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
        outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2L);
        outfile->Printf("\t======================================================================= \n");

        // Get info
        diag = fnoA->occupation_numbers();
        TfnoA = fnoA->nat_orbs();
        VfnoA = fnoA->no_coeff();
        navirA = fnoA->nactvir();
        nfrzv = fnoA->nfrzvir();
        frzvpi_[0] = nfrzv;
        outfile->Printf("\tNumber of active virtuals: %3d\n", navirA);
        outfile->Printf("\tNumber of frozen virtuals: %3d\n", nfrzv);

        // Reset MO coeff
        if (nfrzv > 0) {
            CmoA->copy(VfnoA);
            CoccA.reset();
            CvirA.reset();
            CaoccA.reset();
            CavirA.reset();

            // build mo coeff blocks
            CoccA = std::make_shared<Tensor2d>("Alpha C(mu,i)", nso_, noccA);
            CvirA = std::make_shared<Tensor2d>("Alpha C(mu,a)", nso_, nvirA);
            CaoccA = std::make_shared<Tensor2d>("Alpha Active C(mu,i)", nso_, naoccA);
            CavirA = std::make_shared<Tensor2d>("Alpha Active C(mu,a)", nso_, navirA);
            mo_coeff_blocks();

            // modify pair indices
            navir2AA = navirA * navirA;   // Number of active VIR-VIR pairs
            navoAA = naoccA * navirA;     // Number of active OCC-VIR pairs
            namo = nmo_ - nfrzc - nfrzv;  // Number of active  orbitals
            npop = nmo_ - nfrzv;          // Number of populated orbitals
            ntri_abAA = 0.5 * navirA * (navirA + 1);

            // common things
            common_memfree();
            common_malloc();
        }  // end if nfrzv > 0

        // delete
        fnoA.reset();
        TfnoA.reset();
        VfnoA.reset();
        bQiaA.reset();

    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
    }  // else if (reference_ == "UNRESTRICTED")

    // compute pair indices
    pair_index();

    timer_off("fno");
    outfile->Printf("\tFNOs were generated.\n");
}  // fno_wrapper

}  // namespace dfoccwave
}  // namespace psi
