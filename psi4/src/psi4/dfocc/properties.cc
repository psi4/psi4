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
#include "ekt.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/matrix.h"
#include "psi4/physconst.h"

namespace psi {
namespace dfoccwave {

void DFOCC::oeprop() {
    outfile->Printf("\tComputing one-electron properties...\n");

    timer_on("oeprop");
    // This density finagling won't be necessary if you just set the density on the wavefunction.
    // We do need to set density explicitly on CD-OMP2/OLCCD, at least, since one-electron properties
    // accessible for CD but CD gradients (which would set the density automatically) are not. (That
    // is, CD gradient tech available in the literature but none in Psi4, even for CD-HF.)
    // This should wait until after the current dfocc push. - JPM 08/2022

    auto Da_ = std::make_shared<Matrix>("MO-basis alpha OPDM", nmo_, nmo_);
    auto Db_ = std::make_shared<Matrix>("MO-basis beta OPDM", nmo_, nmo_);
    if (reference_ == "RESTRICTED") {
        G1->to_shared_matrix(Da_);
        Da_->scale(0.5);
        Db_->copy(Da_);
    }

    else if (reference_ == "UNRESTRICTED") {
        G1A->to_shared_matrix(Da_);
        G1B->to_shared_matrix(Db_);
    }

    // Compute oeprop
    auto oe = std::make_shared<OEProp>(shared_from_this());
    oe->set_Da_mo(Da_);
    if (reference_ == "UNRESTRICTED") oe->set_Db_mo(Db_);
    oe->add("DIPOLE");
    oe->add("QUADRUPOLE");
    oe->add("MULLIKEN_CHARGES");
    oe->add("NO_OCCUPATIONS");
    oe->set_title(wfn_type_);
    oe->compute();
    Da_.reset();
    Db_.reset();

    timer_off("oeprop");
}  // end oeprop

//======================================================================
//    EKT-IP
//======================================================================
void DFOCC::ekt_ip() {
    outfile->Printf("\tComputing EKT IPs...\n");

    SharedTensor2d G;
    SharedTensor1d eigA, eigB, psA, psB;

    timer_on("ekt");
    if (reference_ == "RESTRICTED") {

        // Call EKT
        auto ektA = std::make_shared<Ektip>("Alpha EKT", noccA, nmo_, GF, G1, 1.0, 0.5);

        // Print IPs
        outfile->Printf("\n\tEKT Ionization Potentials (Alpha Spin Case) \n");
        outfile->Printf("\t------------------------------------------------------------------- \n");

        // print alpha IPs
        if (print_ < 2) {
            // get only occupieds
            eigA = ektA->eocc();
            psA = ektA->psocc();

            outfile->Printf("\tState    -IP (a.u.)       IP (eV)        Pole Strength \n");
            outfile->Printf("\t------------------------------------------------------------------- \n");
            for (int i = 0; i < noccA; ++i) {
                outfile->Printf("\t%3d %15.6f %15.6f %15.6f \n", i + 1, eigA->get(i), -eigA->get(i) * pc_hartree2ev,
                                psA->get(i));
            }
            outfile->Printf("\t------------------------------------------------------------------- \n");
        }  // end if

        else if (print_ >= 2) {
            // get all
            eigA = ektA->eorb();
            psA = ektA->ps();

            outfile->Printf("\tState    Symmetry   -IP (a.u.)       IP (eV)        Pole Strength \n");
            outfile->Printf("\t------------------------------------------------------------------- \n");
            for (int i = 0; i < noccA; ++i) {
                outfile->Printf("\t%3d %15.6f %15.6f %15.6f \n", i + 1, eigA->get(i), -eigA->get(i) * pc_hartree2ev,
                                psA->get(i));
            }
            outfile->Printf("\t------------------------------------------------------------------- \n");
        }  // end else if

        // delete
        ektA.reset();
        eigA.reset();
        psA.reset();

    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        // Not implemented yet
    }  // else if (reference_ == "UNRESTRICTED")
    timer_off("ekt");
    // outfile->Printf("\tekt is done.\n");
}  // properties.cc

}  // namespace dfoccwave
}  // namespace psi
