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

#include "basewavefunction.h"
#include "wavefunction.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"

namespace psi {

void BaseWavefunction::common_init() {}

BaseWavefunction::BaseWavefunction()
    : options_(Process::environment.options), dipole_field_strength_{{0.0, 0.0, 0.0}}, PCM_enabled_(false) {}

BaseWavefunction::BaseWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis,
                                   Options& options)
    : molecule_(molecule),
      basisset_(basis),
      options_(options),
      dipole_field_strength_{{0.0, 0.0, 0.0}},
      PCM_enabled_(false) {}

BaseWavefunction::BaseWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis)
    : molecule_(molecule),
      basisset_(basis),
      options_(Process::environment.options),
      dipole_field_strength_{{0.0, 0.0, 0.0}},
      PCM_enabled_(false) {}

BaseWavefunction::BaseWavefunction(Options& options)
    : options_(options), dipole_field_strength_{{0.0, 0.0, 0.0}}, PCM_enabled_(false) {}

BaseWavefunction::~BaseWavefunction() = default;

}  // namespace psi
