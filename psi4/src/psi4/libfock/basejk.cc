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

#include "psi4/libfock/basejk.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

namespace psi {

BaseJK::BaseJK(std::shared_ptr<BasisSet> primary) : primary_(primary) { init_knobs(); }

BaseJK::~BaseJK() {}

void BaseJK::init_knobs() {
    print_ = 1;
    debug_ = 0;
    bench_ = 0;

    // 256 MB default
    memory_ = 32000000L;
    omp_nthread_ = 1;
#ifdef _OPENMP
    omp_nthread_ = Process::environment.get_n_threads();
#endif
    cutoff_ = 1.0E-12;

    do_J_ = true;
    do_K_ = true;
    lr_symmetric_ = false;

    num_computed_shells_ = 0L;
    computed_shells_per_iter_ = {};
}

void BaseJK::common_init() {}

const std::unordered_map<std::string, std::vector<size_t>>& BaseJK::computed_shells_per_iter() {
    return computed_shells_per_iter_;
}

const std::vector<size_t>& BaseJK::computed_shells_per_iter(const std::string& n_let) {
    return computed_shells_per_iter_[n_let];
}

size_t BaseJK::num_computed_shells() {
    outfile->Printf(
        "WARNING: BaseJK::num_computed_shells() was called, but benchmarking is disabled for the chosen JK algorithm.");
    outfile->Printf(" Returning 0 as computed shells count.\n");

    return 0;
}

void BaseJK::initialize() { preiterations(); }

void BaseJK::finalize() { postiterations(); }

}  // namespace psi
