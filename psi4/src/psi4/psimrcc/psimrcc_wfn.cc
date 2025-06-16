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

#include "psi4/psifiles.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsio/psio.hpp"

#include "blas.h"
#include "idmrpt2.h"
#include "manybody.h"
#include "mrcc.h"
#include "psimrcc_wfn.h"

namespace psi {

class Options;

namespace psimrcc {
PSIMRCCWfn::PSIMRCCWfn(SharedWavefunction ref_wfn, Options &options) : Wavefunction(options) {
    reference_wavefunction_ = ref_wfn;
    shallow_copy(ref_wfn);
    moinfo_ = std::make_shared<MOInfo>(*(ref_wfn.get()), options);
    moinfo_->setup_model_space();
    memory_ = Process::environment.get_memory();
    free_memory_ = memory_;
}

void PSIMRCCWfn::active_space_warning() const {
    auto nactmo = moinfo_->get_nactv();
    auto nactel = moinfo_->get_nactive_ael() + moinfo_->get_nactive_bel();

    if (nactel > 2 && nactmo > 2) {
        outfile->Printf("\n   WARNING: PSIMRCC detected that you are not using a CAS(2,n) or CAS(m,2) active space");
        outfile->Printf("\n            You requested a CAS(%d,%d) space.  In this case the program will run", nactel,
                        nactmo);
        outfile->Printf("\n            but will negled matrix elements of the effective Hamiltonian between");
        outfile->Printf("\n            reference determinats that differ by more than two spin orbitals.");
        outfile->Printf(
            "\n            The final answer will NOT be the Mk-MRCC energy but only an approximation to it.");
        outfile->Printf("\n            If you are going to report this number in a publication make sure that you");
        outfile->Printf("\n            understand what is going on and that you document it in your publication.");
    }
}

double PSIMRCCWfn::compute_energy() {
    _default_psio_lib_->open(PSIF_PSIMRCC_INTEGRALS, PSIO_OPEN_NEW);

    // Some historical context to understand what you're looking at. PSIMRCC was written around 2006.
    // At the time, it was on the cutting edge of multireference methods. By 2009 and 2010, the main
    // theory-developer decided that Mk-MRCC (which this module exists to compute) wasn't likely to be
    // a viable method, so he stopped development. 2009 was also when the Psi3 -> Psi4 transition started.
    // Accordingly, some aspects of the original module design don't make much sense from the perspective
    // of Psi4, and they were never modernized thoroughly.
    // What you have now is a clumsy attempt to modernize the module from one of the developers brave(?)/foolish
    // enough to hazard dealing with this code. Modernizing here means "uses Wavefunction."
    // Let this serve as a warning against contributing a new module to Psi. New methods should be
    // plugins so you don't create a maintainability problem when you're no longer working on the code.
    // Look on this module, ye mighty, and despair.

    // Ideally, this would go in the constructor, but shared_from_this won't work until after construction.
    if (blas_ == nullptr) {
        blas_ = std::make_shared<CCBLAS>(std::dynamic_pointer_cast<PSIMRCCWfn>(shared_from_this()), options_);
    }

    auto global_timer = std::make_shared<Timer>();
    active_space_warning();

    // TODO: CCManyBody is ancient. Nowadays, they should be wavefunction subclasses.
    std::shared_ptr<CCManyBody> ccmanybody;

    // The astute will notice another method, MP2_CCSD, that is not included here.
    // That is intentional. The method is unpublished, but Francesco is considering finishing it.
    if (options_.get_str("CORR_WFN") == "PT2") {
        ccmanybody = std::make_shared<IDMRPT2>(std::dynamic_pointer_cast<PSIMRCCWfn>(shared_from_this()), options_);
    } else {
        ccmanybody = std::make_shared<CCMRCC>(std::dynamic_pointer_cast<PSIMRCCWfn>(shared_from_this()), options_);
    }

    options_.print();
    auto energy = ccmanybody->compute_energy();

    if (options_.get_str("CORR_WFN") != "PT2") {
        active_space_warning();
    }

    outfile->Printf("\n\n  PSIMRCC job completed.");
    outfile->Printf("\n  Wall Time = %20.6f s", global_timer->get());
    outfile->Printf("\n  GEMM Time = %20.6f s", moinfo_->get_dgemm_timing());

    _default_psio_lib_->close(PSIF_PSIMRCC_INTEGRALS, 1);

    return energy;
}
}
}
