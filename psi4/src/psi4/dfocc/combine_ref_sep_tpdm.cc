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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::combine_ref_sep_tpdm() {
    SharedTensor2d G, Gref, Gsep;
    timer_on("combine_ref_sep_tpdm");
    if (reference_ == "RESTRICTED") {
        //=========================
        // General TPDM
        //=========================
        G = std::make_shared<Tensor2d>("RefSep 3-Index TPDM (Q|AA)", nQ_ref, nmo_, nmo_);

        //=========================
        // OO Block: REF
        //=========================
        Gref = std::make_shared<Tensor2d>("Reference 3-Index TPDM <Q|OO>", nQ_ref, noccA, noccA);
        Gref->read(psio_, PSIF_DFOCC_DENS);
        G->set3_oo(Gref);
        // Gref->print();
        Gref.reset();

        //=========================
        // OO Block: SEP
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM <Q|OO>", nQ_ref, noccA, noccA);
        Gsep->read(psio_, PSIF_DFOCC_DENS);
        G->add3_oo(Gsep, 1.0, 1.0);
        Gsep.reset();

        //=========================
        // OV Block
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM <Q|OV>", nQ_ref, noccA, nvirA);
        Gsep->read(psio_, PSIF_DFOCC_DENS);
        G->set3_ov(Gsep);
        Gsep.reset();

        //=========================
        // VO Block
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM <Q|VO>", nQ_ref, nvirA, noccA);
        Gsep->read(psio_, PSIF_DFOCC_DENS);
        G->set3_vo(Gsep);
        Gsep.reset();

        //=========================
        // VV Block
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM <Q|VV>", nQ_ref, nvirA, nvirA);
        Gsep->read(psio_, PSIF_DFOCC_DENS, true, true);
        G->set3_vv(Gsep, noccA);
        Gsep.reset();

        //=========================
        // Write
        //=========================
        G->write(psio_, PSIF_DFOCC_DENS);
        // G->print();
        G.reset();

    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        //=========================
        // General TPDM: Alpha
        //=========================
        G = std::make_shared<Tensor2d>("RefSep 3-Index TPDM (Q|AA)", nQ_ref, nmo_ * nmo_);

        //=========================
        // OO Block: REF
        //=========================
        Gref = std::make_shared<Tensor2d>("Reference 3-Index TPDM (Q|OO)", nQ_ref, noccA * noccA);
        Gref->read(psio_, PSIF_DFOCC_DENS);
        G->set3_oo(Gref);
        Gref.reset();

        //=========================
        // OO Block: SEP
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|OO)", nQ_ref, noccA * noccA);
        Gsep->read(psio_, PSIF_DFOCC_DENS);
        G->add3_oo(Gsep, 1.0, 1.0);
        Gsep.reset();

        //=========================
        // OV Block
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|OV)", nQ_ref, noccA, nvirA);
        Gsep->read(psio_, PSIF_DFOCC_DENS);
        G->set3_ov(Gsep);
        Gsep.reset();

        //=========================
        // VO Block
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|VO)", nQ_ref, nvirA, noccA);
        Gsep->read(psio_, PSIF_DFOCC_DENS);
        G->set3_vo(Gsep);
        Gsep.reset();

        //=========================
        // VV Block
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|VV)", nQ_ref, nvirA, nvirA);
        Gsep->read(psio_, PSIF_DFOCC_DENS, true, true);
        G->set3_vv(Gsep, noccA);
        Gsep.reset();

        //=========================
        // Write
        //=========================
        G->write(psio_, PSIF_DFOCC_DENS);
        G->print();
        G.reset();

        //=========================
        // General TPDM: Beta
        //=========================
        G = std::make_shared<Tensor2d>("RefSep 3-Index TPDM (Q|aa)", nQ_ref, nmo_ * nmo_);

        //=========================
        // oo Block: REF
        //=========================
        Gref = std::make_shared<Tensor2d>("Reference 3-Index TPDM (Q|oo)", nQ_ref, noccB * noccB);
        Gref->read(psio_, PSIF_DFOCC_DENS);
        G->set3_oo(Gref);
        Gref.reset();

        //=========================
        // oo Block: SEP
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|oo)", nQ_ref, noccB * noccB);
        Gsep->read(psio_, PSIF_DFOCC_DENS);
        G->add3_oo(Gsep, 1.0, 1.0);
        Gsep.reset();

        //=========================
        // ov Block
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|ov)", nQ_ref, noccB, nvirB);
        Gsep->read(psio_, PSIF_DFOCC_DENS);
        G->set3_ov(Gsep);
        Gsep.reset();

        //=========================
        // vo Block
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|vo)", nQ_ref, nvirB, noccB);
        Gsep->read(psio_, PSIF_DFOCC_DENS);
        G->set3_vo(Gsep);
        Gsep.reset();

        //=========================
        // vv Block
        //=========================
        Gsep = std::make_shared<Tensor2d>("3-Index Separable TPDM (Q|vv)", nQ_ref, nvirB, nvirB);
        Gsep->read(psio_, PSIF_DFOCC_DENS, true, true);
        G->set3_vv(Gsep, noccB);
        Gsep.reset();

        //=========================
        // Write
        //=========================
        G->write(psio_, PSIF_DFOCC_DENS);
        G->print();
        G.reset();

    }  // else if (reference_ == "UNRESTRICTED")
    timer_off("combine_ref_sep_tpdm");
}  // end combine_ref_sep_tpdm

}  // namespace dfoccwave
}  // namespace psi
