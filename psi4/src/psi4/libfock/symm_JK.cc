/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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
#include "psi4/libmints/sieve.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/lib3index/df_helper.h"

#include "jk.h"

#include <sstream>
#include "psi4/libpsi4util/PsiOutStream.h"
#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {

symm_JK::symm_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary)
    : JK(primary), auxiliary_(auxiliary) {
    common_init();
}

symm_JK::~symm_JK() {}

void symm_JK::common_init() {
    
    dfh_ = std::make_shared<DF_Helper>(primary_, auxiliary_);

}

void symm_JK::preiterations() {
   
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
 
    // we need to prepare the AOs here, and that's it.
    // DF_Helper takes care of all the housekeeping

    // TODO add wK integrals.
    if (!do_wK_) {
        dfh_->initialize();
    } else {
           ; // initialize_wK()  
    }

}
void symm_JK::compute_JK() {

    dfh_->build_JK(C_left_ao_, C_right_ao_, D_ao_, J_ao_, K_ao_, 
                   max_nocc(), do_J_, do_K_, do_wK_);

}
void symm_JK::postiterations() {
}
void symm_JK::print_header() const {
    dfh_->print_header();
}
int symm_JK::max_nocc() const {
    int max_nocc = 0;
    for (size_t N = 0; N < C_left_ao_.size(); N++) {
        max_nocc = (C_left_ao_[N]->colspi()[0] > max_nocc ? C_left_ao_[N]->colspi()[0] : max_nocc);
    }
    return max_nocc;
}
}
