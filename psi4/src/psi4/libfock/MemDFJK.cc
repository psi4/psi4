/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/lib3index/dfhelper.h"

#include "jk.h"

#include <sstream>
#include "psi4/libpsi4util/PsiOutStream.h"
#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {

MemDFJK::MemDFJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options) : JK(primary), auxiliary_(auxiliary), options_(options) {
    common_init();
}

MemDFJK::~MemDFJK() {}

void MemDFJK::common_init() { dfh_ = std::make_shared<DFHelper>(primary_, auxiliary_, options_); }
size_t MemDFJK::memory_estimate() {
    dfh_->set_nthreads(omp_nthread_);
    dfh_->set_schwarz_cutoff(cutoff_);
    return dfh_->get_core_size();
}

void MemDFJK::preiterations() {
    // Initialize calls your derived class's preiterations member
    // knobs are set and state variables assigned

    // use previously set state variables to dfh instance
    dfh_->set_nthreads(omp_nthread_);
    dfh_->set_schwarz_cutoff(cutoff_);
    dfh_->set_method("STORE");
    dfh_->set_fitting_condition(condition_);
    dfh_->set_memory(memory_ - memory_overhead());
    dfh_->set_do_wK(do_wK_);
    dfh_->set_omega(omega_);
    if (do_wK_) { 
        dfh_->set_wcombine(wcombine_); 
    } else {
        dfh_->set_wcombine(false);
        wcombine_ = false;
    }
    dfh_->set_omega_alpha(omega_alpha_);
    dfh_->set_omega_beta(omega_beta_);

    // we need to prepare the AOs here, and that's it.
    // DFHelper takes care of all the housekeeping

    dfh_->initialize();
}
void MemDFJK::compute_JK() {

    // zero out J, K, and wK matrices
    zero();

    dfh_->build_JK(C_left_ao_, C_right_ao_, D_ao_, J_ao_, K_ao_, wK_ao_, max_nocc(), do_J_, do_K_, do_wK_,
                   lr_symmetric_);
    if (lr_symmetric_) {
        if (do_wK_) {
            for (size_t N = 0; N < wK_ao_.size(); N++) {
                wK_ao_[N]->hermitivitize();
            }
        }
    }
}
void MemDFJK::postiterations() {}
void MemDFJK::print_header() const {
    // dfh_->print_header();
    if (print_) {
        outfile->Printf("  ==> MemDFJK: Density-Fitted J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:           %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:           %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:          %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:              %11.3E\n", omega_);
        outfile->Printf("    OpenMP threads:     %11d\n", omp_nthread_);
        outfile->Printf("    Memory [MiB]:       %11ld\n", (memory_ * 8L) / (1024L * 1024L));
        outfile->Printf("    Algorithm:          %11s\n", (dfh_->get_AO_core() ? "Core" : "Disk"));
        outfile->Printf("    Schwarz Cutoff:     %11.0E\n", cutoff_);
        outfile->Printf("    Mask sparsity (%%):  %11.4f\n", 100. * dfh_->ao_sparsity());
        outfile->Printf("    Fitting Condition:  %11.0E\n\n", condition_);

        outfile->Printf("   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }
}
int MemDFJK::max_nocc() const {
    int max_nocc = 0;
    for (size_t N = 0; N < C_left_ao_.size(); N++) {
        max_nocc = (C_left_ao_[N]->colspi()[0] > max_nocc ? C_left_ao_[N]->colspi()[0] : max_nocc);
    }
    return max_nocc;
}
void MemDFJK::set_omega_alpha(double alpha) {
    omega_alpha_ = alpha;
    dfh_->set_omega_alpha(omega_alpha_);
}
void MemDFJK::set_omega_beta(double beta){
    omega_beta_ = beta;
    dfh_->set_omega_beta(omega_beta_);
}
void MemDFJK::set_do_wK(bool tf) { do_wK_ = tf; dfh_->set_do_wK(tf); }
void MemDFJK::set_wcombine(bool wcombine) { 
    wcombine_ = wcombine;
    if (dfh_) {
        dfh_->set_wcombine(wcombine); 
    }
}
}  // namespace psi
