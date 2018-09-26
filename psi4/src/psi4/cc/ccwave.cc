/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include "ccwave.h"

#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

namespace psi {
namespace cc {

CCWavefunction::CCWavefunction(std::shared_ptr<Wavefunction> ref_wfn) : Wavefunction(Process::environment.options) {
    // Copy the wavefuntion then update
    shallow_copy(ref_wfn);
    set_reference_wavefunction(ref_wfn);
    common_init();
}

CCWavefunction::CCWavefunction(std::shared_ptr<Wavefunction> ref_wfn, Options &options) : Wavefunction(options) {
    // Copy the wavefuntion then update
    shallow_copy(ref_wfn);
    set_reference_wavefunction(ref_wfn);
    common_init();
}

CCWavefunction::~CCWavefunction() {}

void CCWavefunction::common_init() {}

double CCWavefunction::compute_energy() { return 0.0; }

}  // namespace cc
}  // namespace psi
