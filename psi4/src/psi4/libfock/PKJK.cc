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

#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "jk.h"
#include "PKmanagers.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/sobasis.h"

#include <sstream>
#include "psi4/libpsi4util/PsiOutStream.h"
#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {

PKJK::PKJK(std::shared_ptr<BasisSet> primary, Options& options) : JK(primary), options_(options) { common_init(); }

PKJK::~PKJK() {}

void PKJK::common_init() {
    pk_file_ = PSIF_SO_PK;
    nthreads_ = 1;
#ifdef _OPENMP
    nthreads_ = Process::environment.get_n_threads();
#endif
}
size_t PKJK::memory_estimate() {
    size_t nbf = (size_t)primary_->nbf();
    return nbf * nbf * nbf * nbf;
}

bool PKJK::C1() const { return true; }

void PKJK::print_header() const {
    if (print_) {
        outfile->Printf("  ==> DiskJK: Disk-Based J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    Memory [MiB]:      %11ld\n", (memory_ * 8L) / (1024L * 1024L));
        outfile->Printf("    Schwarz Cutoff:    %11.0E\n\n", cutoff_);
        outfile->Printf("    OpenMP threads:    %11d\n\n", nthreads_);
    }
}

void PKJK::preiterations() {
    // Build PKManager to get proper algorithm set up
    Options& options = Process::environment.options;

    psio_ = _default_psio_lib_;

    timer_on("Total PK formation time");
    // We compute the integrals so that we can directly write the
    // PK file to disk. Also, do everything in the AO basis
    // like the modern JK algos, for adding sieving later

    PKmanager_ = pk::PKManager::build_PKManager(psio_, primary_, memory_, options, do_wK_, omega_);

    PKmanager_->initialize();

    PKmanager_->form_PK();

    // If range-separated K needed, we redo all the above steps
    if (do_wK_) {
        outfile->Printf("  Computing range-separated integrals for PK\n");

        PKmanager_->initialize_wK();

        PKmanager_->form_PK_wK();
    }

    // PK files are written at this point. We are done.
    timer_off("Total PK formation time");
}

void PKJK::compute_JK() {

    // zero out J, K, and wK matrices
    zero();

    timer_on("PK computes JK");
    // We form the vector containing the density matrix triangular elements
    PKmanager_->prepare_JK(D_ao_, C_left_ao_, C_right_ao_);

    if (J_ao_.size()) {
        // We can safely pass K here since its size is checked within
        // the routine
        PKmanager_->form_J(J_ao_, "", K_ao_);
    }
    if (K_ao_.size()) {
        PKmanager_->form_K(K_ao_);
    }
    if (wK_ao_.size()) {
        PKmanager_->form_wK(wK_ao_);
    }

    PKmanager_->finalize_JK();

    timer_off("PK computes JK");
}

void PKJK::postiterations() {}
}  // namespace psi
