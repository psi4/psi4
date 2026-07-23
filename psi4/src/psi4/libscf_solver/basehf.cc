/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "basehf.h"

#include "psi4/libmints/molecule.h"
#include "psi4/liboptions/liboptions.h"

namespace psi {
namespace scf {

void BaseHF::common_init(Options& options, std::string& module, const std::shared_ptr<Molecule>& molecule,
                         const std::array<double, 3>& dipole_field_strength) {
    attempt_number_ = 1;
    reset_occ_ = false;
    sad_ = false;
    module = "scf";

    scf_type_ = options.get_str("SCF_TYPE");

    energies_["Total Energy"] = 0.0;
    energies_["-D"] = 0.0;

    nuclearrep_ = molecule->nuclear_repulsion_energy(dipole_field_strength);
    charge_ = molecule->molecular_charge();
    multiplicity_ = molecule->multiplicity();

    initialized_diis_manager_ = false;
    MOM_performed_ = false;  // duplicated py-side (needed before iterate)
}

}  // namespace scf
}  // namespace psi
