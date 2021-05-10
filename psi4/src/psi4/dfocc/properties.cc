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
#include "ekt.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/dipole.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/physconst.h"

namespace psi {
namespace dfoccwave {

void DFOCC::oeprop() {
    outfile->Printf("\tComputing one-electron properties...\n");

    timer_on("oeprop");
    // This density finagling won't be necessary if you just set the density on the wavefunction.
    // However, CD-OMP2/OLCCD have oeprop enabled but gradients disabled for undocumented reasons.
    // Accordingly, we'll play it safe and not set the density for those cases just yet.
    // ...Which means we can't assume the density to be set in general. Sad. - JPM 07/2020
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

    // Compute dipole moments
    Vector3 origin;
    origin[0] = 0.0;
    origin[1] = 0.0;
    origin[2] = 0.0;
    SharedVector ndip = DipoleInt::nuclear_contribution(molecule_, origin);

    // Compute dipoles
    std::vector<SharedMatrix> Xtemp;
    SharedTensor2d Xso_, Yso_, Zso_, Gtemp;
    Xso_ = std::make_shared<Tensor2d>("Dipole X <mu|nu>", nso_, nso_);
    Yso_ = std::make_shared<Tensor2d>("Dipole Y <mu|nu>", nso_, nso_);
    Zso_ = std::make_shared<Tensor2d>("Dipole Z <mu|nu>", nso_, nso_);
    Gtemp = std::make_shared<Tensor2d>("G AO <mu|nu>", nso_, nso_);

    // Read SO-basis one-electron integrals
    MintsHelper mints(MintsHelper(reference_wavefunction_->basisset(), options_, 0));
    Xtemp = mints.ao_dipole();
    for(int mu = 0; mu < nso_; ++mu) {
        for(int nu = 0; nu < nso_; ++nu) {
            Xso_->set(mu,nu, Xtemp[0]->get(mu,nu));
            Yso_->set(mu,nu, Xtemp[1]->get(mu,nu));
            Zso_->set(mu,nu, Xtemp[2]->get(mu,nu));
        }
    }
    //Xso_->print();

    // Compute dipole moments
    if (reference_ == "RESTRICTED") {
        Gtemp->back_transform(G1, CmoA);
    }       
    else if (reference_ == "UNRESTRICTED") {
        Gtemp->back_transform(G1A, CmoA);
        Gtemp->back_transform(G1B, CmoB, 1.0, 1.0);
    }
    double Mux = ndip->get(0) + Gtemp->vector_dot(Xso_); 
    double Muy = ndip->get(1) + Gtemp->vector_dot(Yso_); 
    double Muz = ndip->get(2) + Gtemp->vector_dot(Zso_); 

    // Open file for writing.
    outfile->Printf("\tMu_X (a.u.): %12.12f\n", Mux);
    outfile->Printf("\tMu_Y (a.u.): %12.12f\n", Muy);
    outfile->Printf("\tMu_Z (a.u.): %12.12f\n", Muz);
    Xso_.reset();
    Yso_.reset();
    Zso_.reset();
    Gtemp.reset();

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
        // malloc
        eigA = std::make_shared<Tensor1d>("epsilon <I|J>", noccA);
        psA = std::make_shared<Tensor1d>("alpha occupied pole strength vector", noccA);

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
