/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
#include "psi4/libciomr/libciomr.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libqt/qt.h"

#include "blas.h"
#include "main.h"
#include "sort.h"
#include "mrcc.h"
#include "psimrcc_wfn.h"
#include "transform.h"

// PSI FILES

using namespace psi;

namespace psi {

namespace psimrcc {
// Global variables - created in compute_energy

SharedWavefunction psimrcc(SharedWavefunction ref_wfn, Options &options) {
    using namespace psi::psimrcc;

    outfile->Printf("\n  MRCC          MRCC");
    outfile->Printf("\n   MRCC  MRCC  MRCC");
    outfile->Printf("\n   MRCC  MRCC  MRCC      PSIMRCC Version 0.9.3.3, July 2009");
    outfile->Printf("\n   MRCC  MRCC  MRCC      Multireference Coupled Cluster, written by");
    outfile->Printf("\n     MRCCMRCCMRCC        Francesco A. Evangelista and Andrew C. Simmonett");
    outfile->Printf("\n         MRCC            Compiled on %s at %s", __DATE__, __TIME__);
    outfile->Printf("\n         MRCC");
    outfile->Printf("\n       MRCCMRCC");

    auto wfn = std::make_shared<PSIMRCCWfn>(ref_wfn, options);

    if (options["PERTURB_CBS"].has_changed() || options["PERTURB_CBS_COUPLING"].has_changed()) {
        outfile->Printf("\tPerturbative CBS was removed in 1.4. Using unpublished features is a bad habit.\n\n");
    }

    wfn->compute_energy();

    return wfn;
}

}  // namespace psimrcc
}  // namespace psi
