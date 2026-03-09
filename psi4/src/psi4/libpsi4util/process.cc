/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include "psi4/psi4-dec.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/extern.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

// MKL Header
#ifdef USING_LAPACK_MKL
#include <mkl.h>
#endif

// OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace psi {

Process::Environment Process::environment;
const std::string empty_;

void Process::Environment::initialize() {
    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

    /* See notes in process.h */
    _psio_manager_keepalive = PSIOManager::shared_object();
}

void Process::Environment::set_n_threads(int nthread) {
    nthread_ = nthread;
#ifdef _OPENMP
    omp_set_num_threads(nthread_);
#endif
#ifdef USING_LAPACK_MKL
    mkl_set_num_threads(nthread_);
#endif

    // HACK: TODO: CC-pthread codes should ask us how many threads
    // No, this didn't work, back this out for now (and we won't need
    // this once the final solution is in anyway)  --CDS
    // Process::environment.options.set_global_int("NUM_THREADS",nthread);
}

void Process::Environment::set_molecule(const std::shared_ptr<Molecule> &molecule) { molecule_ = molecule; }

std::shared_ptr<Molecule> Process::Environment::molecule() const { return molecule_; }

Process::Environment Process::get_environment() { return environment; }

size_t Process::Environment::get_memory() const { return memory_; }

void Process::Environment::set_memory(size_t m) { memory_ = m; }

int Process::Environment::get_n_threads() const { return nthread_; }

void die_if_not_converged() {
    outfile->Printf("Iterations did not converge.");

    if (Process::environment.options.get_bool("DIE_IF_NOT_CONVERGED"))
        throw PSIEXCEPTION("Iterations did not converge.");
    else {
        outfile->Printf("Iterations did not converge.");
    }
}
}
