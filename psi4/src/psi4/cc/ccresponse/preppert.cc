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

/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

void write_blocks(const Matrix& mat);

/* preppert(): Prepare DPD structures for all currently needed one-electron
** property integrals in the MO basis.
**
** -TDC, 6/11
*/

void preppert(std::shared_ptr<BasisSet> primary) {
    std::vector<std::string> cartstr {"X", "Y", "Z"};

    MintsHelper mints(primary, Process::environment.options, 0);
    auto dipole = mints.so_dipole();
    auto nabla = mints.so_nabla();
    auto angmom = mints.so_angular_momentum();
    auto trquad = mints.so_traceless_quadrupole();

    // Electric dipole integrals
    for (int i = 0; i < 3; i++) {
        dipole[i]->transform(moinfo.Ca);
        dipole[i]->set_name("Mu_" + cartstr[i]);
        write_blocks(*dipole[i]);
    }

    // Velocity-gauge electric dipole integrals
    for (int i = 0; i < 3; i++) {
        nabla[i]->transform(moinfo.Ca);
        nabla[i]->set_name("P_" + cartstr[i]);
        write_blocks(*nabla[i]);
    }

    // Complex conjugate of velocity-gauge electric dipole integrals
    for (int i = 0; i < 3; i++) {
        nabla[i]->scale(-1.0);
        nabla[i]->set_name("P*_" + cartstr[i]);
        write_blocks(*nabla[i]);
    }

    // Magnetic dipole integrals (these require a -1/2 prefactor)
    for (int i = 0; i < 3; i++) {
        angmom[i]->scale(-0.5);
        angmom[i]->transform(moinfo.Ca);
        angmom[i]->set_name("L_" + cartstr[i]);
        write_blocks(*angmom[i]);
    }

    // Complex conjugate of magnetic dipole integrals
    for (int i = 0; i < 3; i++) {
        angmom[i]->scale(-1.0);
        angmom[i]->set_name("L*_" + cartstr[i]);
        write_blocks(*angmom[i]);
    }

    // Traceless quadrupole integrals
    for (int i = 0, ij = 0; i < 3; i++) {
        for (int j = i; j < 3; j++, ij++) {
            trquad[ij]->transform(moinfo.Ca);
            trquad[ij]->set_name("Q_" + cartstr[i] + cartstr[j]);
            write_blocks(*trquad[ij]);
            if (i != j) {
                trquad[ij]->set_name("Q_" + cartstr[j] + cartstr[i]);
                write_blocks(*trquad[ij]);
            }
        }
    }
}

}  // namespace ccresponse
}  // namespace psi
