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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"
#include "fno.h"
#include "psi4/libmints/matrix.h"
#include "psi4/physconst.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"

namespace psi {
namespace dfoccwave {

void DFOCC::fno_wrapper() {
    outfile->Printf("\n\tComputing FNOs...\n");

    SharedTensor2d T;
    SharedTensor1d diag, diagA, diagB;

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

    if (reference_ == "UNRESTRICTED") {
        bQiaB = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB);
        bQiaB->read(psio_, PSIF_DFOCC_INTS);
    }

    // Remove old integrals file
    remove_binary_file(PSIF_DFOCC_INTS);
    // outfile->Printf("\tI am here\n");

    timer_on("fno");
    //===========================//
    //=== RHF FNO ===============//
    //===========================//
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

        Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2L;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Ecorr;

        // Get info
        diag = fnoA->occupation_numbers();
        TfnoA = fnoA->nat_orbs();
        VfnoA = fnoA->no_coeff();
        navirA = fnoA->nactvir();
        nfrzv = fnoA->nfrzvir();
        frzvpi_[0] = nfrzv;
        outfile->Printf("\tNumber of active virtuals: %3d\n", navirA);
        outfile->Printf("\tNumber of frozen virtuals: %3d\n", nfrzv);
        
        //VfnoA->print();

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

    //===========================//
    //=== UHF FNO ===============//
    //===========================//
    else if (reference_ == "UNRESTRICTED") {
        // malloc
        TfnoA = std::make_shared<Tensor2d>("FNO Alpha T matrix", nmo_, nmo_);
        TfnoB = std::make_shared<Tensor2d>("FNO Beta T matrix", nmo_, nmo_);
        VfnoA = std::make_shared<Tensor2d>("FNO Alpha V matrix", nso_, nmo_);
        VfnoB = std::make_shared<Tensor2d>("FNO Beta V matrix", nso_, nmo_);
        diagA = std::make_shared<Tensor1d>("Alpha Occupation numbers", nvirA);
        diagB = std::make_shared<Tensor1d>("Beta Occupation numbers", nvirB);

        // Call FnoCC
        SharedFnoCC fno;
        int navir_fno = 0;
        if (options_["ACTIVE_NAT_ORBS"].has_changed()) {
            for (int h = 0; h < nirrep_; h++) {
                navir_fno += (int)options_["ACTIVE_NAT_ORBS"][h].to_double();
            }
            fno = std::make_shared<FnoCC>("UHF FNO", nQ, nso_, nfrzc, naoccA, naoccB, nvirA, nvirB, 
                                           CmoA, CmoB, FockA, FockB, bQiaA, bQiaB, 
                                           tol_fno, fno_percentage, navir_fno, false, true);
        }

        else if (options_["OCC_PERCENTAGE"].has_changed()) {
            fno = std::make_shared<FnoCC>("UHF FNO", nQ, nso_, nfrzc, naoccA, naoccB, nvirA, nvirB, 
                                           CmoA, CmoB, FockA, FockB, bQiaA, bQiaB, 
                                           tol_fno, fno_percentage, navir_fno, true, false);

        }

        else {
            fno = std::make_shared<FnoCC>("UHF FNO", nQ, nso_, nfrzc, naoccA, naoccB, nvirA, nvirB, 
                                           CmoA, CmoB, FockA, FockB, bQiaA, bQiaB, 
                                           tol_fno, fno_percentage, navir_fno, false, false);
        }

        // Get MP2 energy
        Ecorr = fno->mp2_corr();
        Emp2L = Eref + Ecorr;
        outfile->Printf("\tComputing DF-MP2 energy in the full virtual space... \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tDF-HF Energy (a.u.)                : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tDF-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
        outfile->Printf("\tDF-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2L);
        outfile->Printf("\t======================================================================= \n");

        Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2L;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = Ecorr;

        // Get info
        diagA = fno->occupation_numbers_alpha();
        diagB = fno->occupation_numbers_beta();
        TfnoA = fno->nat_orbs_alpha();
        TfnoB = fno->nat_orbs_beta();
        VfnoA = fno->no_coeff_alpha();
        VfnoB = fno->no_coeff_beta();

        navirA = fno->nactvir_alpha();
        navirB = fno->nactvir_beta();
        nfrzv = fno->nfrzvir();
        int nfrzvA = fno->nfrzvir_alpha();
        int nfrzvB = fno->nfrzvir_beta();
        frzvpi_[0] = nfrzv;
        outfile->Printf("\tNumber of active alpha virtuals: %3d\n", navirA);
        outfile->Printf("\tNumber of active beta virtuals : %3d\n", navirB);
        //outfile->Printf("\tNumber of frozen alpha virtuals: %3d\n", nfrzvA);
        //outfile->Printf("\tNumber of frozen beta virtuals : %3d\n", nfrzvB);
        outfile->Printf("\tNumber of frozen virtuals      : %3d\n", nfrzv);

        //VfnoA->print();
        //VfnoB->print();

        // Reset MO coeff
        if (nfrzv > 0) {
            CmoA->copy(VfnoA);
            CoccA.reset();
            CvirA.reset();
            CaoccA.reset();
            CavirA.reset();

            CmoB->copy(VfnoB);
            CoccB.reset();
            CvirB.reset();
            CaoccB.reset();
            CavirB.reset();

            // build mo coeff blocks
            CoccA = std::make_shared<Tensor2d>("Alpha C(mu,i)", nso_, noccA);
            CvirA = std::make_shared<Tensor2d>("Alpha C(mu,a)", nso_, nvirA);
            CaoccA = std::make_shared<Tensor2d>("Alpha Active C(mu,i)", nso_, naoccA);
            CavirA = std::make_shared<Tensor2d>("Alpha Active C(mu,a)", nso_, navirA);

            CoccB = std::make_shared<Tensor2d>("Beta C(mu,i)", nso_, noccB);
            CvirB = std::make_shared<Tensor2d>("Beta C(mu,a)", nso_, nvirB);
            CaoccB = std::make_shared<Tensor2d>("Beta Active C(mu,i)", nso_, naoccB);
            CavirB = std::make_shared<Tensor2d>("Beta Active C(mu,a)", nso_, navirB);
            mo_coeff_blocks();

            // modify pair indices
            navir2AA = navirA * navirA;   // Number of active VIR-VIR pairs
            navoAA = naoccA * navirA;     // Number of active OCC-VIR pairs
            namo = nmo_ - nfrzc - nfrzv;  // Number of active  orbitals
            npop = nmo_ - nfrzv;          // Number of populated orbitals
            ntri_abAA = navirA * (navirA + 1) / 2;

            navir2BB = navirB * navirB;   // Number of active VIR-VIR pairs
            navoBB = naoccB * navirB;     // Number of active OCC-VIR pairs
            //namo = nmo_ - nfrzc - nfrzv;  // Number of active  orbitals
            //npop = nmo_ - nfrzv;          // Number of populated orbitals
            ntri_abBB = navirB * (navirB + 1) / 2;

            if (naoccA > 1)
                ntri_anti_ijAA = 0.5 * naoccA * (naoccA - 1);
            else
                ntri_anti_ijAA = 1;
            if (naoccB > 1)
                ntri_anti_ijBB = 0.5 * naoccB * (naoccB - 1);
            else
                ntri_anti_ijBB = 1;
            if (navirA > 1)
                ntri_anti_abAA = 0.5 * navirA * (navirA - 1);
            else
                ntri_anti_abAA = 1;
            if (navirB > 1)
                ntri_anti_abBB = 0.5 * navirB * (navirB - 1);
            else
                ntri_anti_abBB = 1;

            // common things
            common_memfree();
            common_malloc();
        }  // end if nfrzv > 0

        // delete
        fno.reset();
        TfnoA.reset();
        TfnoB.reset();
        VfnoA.reset();
        VfnoB.reset();
        bQiaA.reset();
        bQiaB.reset();

    }  // else if (reference_ == "UNRESTRICTED")

    // compute pair indices
    pair_index();

    timer_off("fno");
    outfile->Printf("\tFNOs were generated.\n");
}  // fno_wrapper

}  // namespace dfoccwave
}  // namespace psi
