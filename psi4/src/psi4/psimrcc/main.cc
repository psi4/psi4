/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

/*******************************************************************
 *  PSIMRCC (2008) by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ******************************************************************/

/**
 *  @defgroup PSIMRCC PSIMRCC is a code for SR/MRCC computations
 *  @file psimrcc.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Contains main() and global variables
 */

#include "psi4/psi4-dec.h"
// Standard libraries
#include <iostream>
#include <complex>
#include <cstdlib>

// PSI libraries
#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/memory_manager.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"

// PSI C++
#include "psi4/libpsio/psio.hpp"

#include "blas.h"
#include "main.h"
#include "sort.h"
#include "mrcc.h"
#include "psimrcc.h"
#include "transform.h"

// PSI FILES

using namespace psi;

namespace psi {

namespace psimrcc {
// Global variables
CCBLAS *blas;
CCSort *sorter;
CCTransform *trans = nullptr;
MOInfo *moinfo;
ModelSpace *model_space;
MemoryManager *memory_manager;

PsiReturnType psimrcc(SharedWavefunction ref_wfn, Options &options) {
    using namespace psi::psimrcc;
    tstart();
    _default_psio_lib_->open(PSIF_PSIMRCC_INTEGRALS, PSIO_OPEN_NEW);

    outfile->Printf("\n  MRCC          MRCC");
    outfile->Printf("\n   MRCC  MRCC  MRCC");
    outfile->Printf("\n   MRCC  MRCC  MRCC      PSIMRCC Version 0.9.3.3, July 2009");
    outfile->Printf("\n   MRCC  MRCC  MRCC      Multireference Coupled Cluster, written by");
    outfile->Printf("\n     MRCCMRCCMRCC        Francesco A. Evangelista and Andrew C. Simmonett");
    outfile->Printf("\n         MRCC            Compiled on %s at %s", __DATE__, __TIME__);
    outfile->Printf("\n         MRCC");
    outfile->Printf("\n       MRCCMRCC");

    moinfo = new MOInfo(*(ref_wfn.get()), options);

    memory_manager = new MemoryManager(Process::environment.get_memory());
    model_space = new ModelSpace(moinfo);

    moinfo->setup_model_space();  // The is a bug here DELETEME

    if (options["PERTURB_CBS"].has_changed() || options["PERTURB_CBS_COUPLING"].has_changed()) {
        outfile->Printf("\tPerturbative CBS was removed in 1.4. Using unpublished features is a bad habit.\n\n");
    }

    int nactmo = moinfo->get_nactv();
    int nactel = moinfo->get_nactive_ael() + moinfo->get_nactive_bel();
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

    blas = new CCBLAS(options);
    trans = new CCTransform(ref_wfn);
    if (options.get_str("CORR_WFN") == "PT2") {
        mrpt2(ref_wfn, options);
    } else {
        mrccsd(ref_wfn, options);
        if (nactel > 2 && nactmo > 2) {
            outfile->Printf(
                "\n   WARNING: PSIMRCC detected that you are not using a CAS(2,n) or CAS(m,2) active space");
            outfile->Printf("\n            You requested a CAS(%d,%d) space.  In this case the program will run",
                            nactel, nactmo);
            outfile->Printf("\n            but will negled matrix elements of the effective Hamiltonian between");
            outfile->Printf("\n            reference determinats that differ by more than two spin orbitals.");
            outfile->Printf(
                "\n            The final answer will NOT be the Mk-MRCC energy but only an approximation to it.");
            outfile->Printf("\n            If you are going to report this number in a publication make sure that you");
            outfile->Printf("\n            understand what is going on and that you document it in your publication.");
        }
    }

    delete sorter;
    delete trans;
    delete blas;

    memory_manager->MemCheck("outfile");

    delete model_space;
    delete moinfo;
    delete memory_manager;

    _default_psio_lib_->close(PSIF_PSIMRCC_INTEGRALS, 1);
    tstop();

    return Success;
}

}  // namespace psimrcc
}  // namespace psi
