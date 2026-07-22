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

#include "cghf.h"

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libpsi4util/process.h"

namespace psi {
namespace scf {

CGHF::CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> func)
    : CGHF(ref_wfn, func, Process::environment.options, PSIO::shared_object()) {}

CGHF::CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> func, Options& options,
           std::shared_ptr<PSIO> psio)
    : ComplexWavefunction(options), functional_(func) {
    // Analogous to HF::HF's shallow_copy(ref_wfn); full ComplexWavefunction::shallow_copy comes later.
    molecule_ = ref_wfn->molecule();
    basisset_ = ref_wfn->basisset();
    psio_ = psio;
    // Options-only base ctor skipped init; run it now that molecule/basis are set.
    ComplexWavefunction::common_init();
    common_init();
}

CGHF::~CGHF() {}

void CGHF::common_init() {
    name_ = "CGHF";
    module_ = "scf";
}

double CGHF::compute_energy() {
    // TODO: Implement the CGHF SCF energy (initialize / iterate / finalize).
    return 0.0;
}

}  // namespace scf
}  // namespace psi
