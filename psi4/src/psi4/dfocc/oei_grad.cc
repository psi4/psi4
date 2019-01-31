/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
#include <fstream>
#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/mintshelper.h"
#include "dfocc.h"

#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::oei_grad() {
    /********************************************************************************************/
    /************************** Gradient ********************************************************/
    /********************************************************************************************/
    // outfile->Printf("\tComputing analytic gradients...\n");
    //

    /********************************************************************************************/
    /************************** Nuclear Gradient ************************************************/
    /********************************************************************************************/
    // => Nuclear Gradient <= //
    gradients["Nuclear"] = SharedMatrix(molecule_->nuclear_repulsion_energy_deriv1(dipole_field_strength_).clone());
    gradients["Nuclear"]->set_name("Nuclear Gradient");
    gradients["Nuclear"]->print_atom_vector();

    // => Kinetic Gradient <= //
    auto mints = std::make_shared<MintsHelper>(basisset_, options_);

    /****************************************************************************************/
    /************************** One-Electron Gradient ***************************************/
    /****************************************************************************************/
    timer_on("Grad: V T Perturb");
    auto D = std::make_shared<Matrix>("AO-basis OPDM", nmo_, nmo_);
    G1ao->to_shared_matrix(D);
    gradients["Core"] = mints->core_hamiltonian_grad(D);
    gradients["Core"]->print_atom_vector();
    timer_off("Grad: V T Perturb");

    /****************************************************************************************/
    /************************** Overlap Gradient ********************************************/
    /****************************************************************************************/
    // => Overlap Gradient <= //
    timer_on("Grad: S");
    auto W = std::make_shared<Matrix>("AO-basis Energy-Weighted OPDM", nmo_, nmo_);
    GFao->to_shared_matrix(W);
    gradients["Overlap"] = mints->overlap_grad(W);
    gradients["Overlap"]->scale(-1.0);
    gradients["Overlap"]->print_atom_vector();
    timer_off("Grad: S");

    /****************************************************************************************/
    /****************************************************************************************/
    /****************************************************************************************/

    // outfile->Printf("\toei_grad is done. \n");
}
}  // namespace dfoccwave
}  // namespace psi
