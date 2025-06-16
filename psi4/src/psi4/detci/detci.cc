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

/*! \defgroup DETCI detci: The Determinant CI code */

/*! \file
    \ingroup DETCI
    \brief Determinant-based CI program

   DETCI

   DETERMINANT CI Program, incorporating Abelian point-group symmetry

   C. David Sherrill
   Center for Computational Quantum Chemistry
   University of Georgia
   August 1994

   Updated 3/95 to do frozen core and virtuals correctly
   Updated 5/95 to do RAS CI's again
   Updated 2/96 to clean up code and rename DETCI

*/

#include <cstdio>

#include "psi4/detci/structs.h"
#include "psi4/detci/ciwave.h"

namespace psi {
namespace detci {

SharedWavefunction detci(SharedWavefunction ref_wfn, Options &options);
}
}  // namespace psi

namespace psi {
namespace detci {

SharedWavefunction detci(SharedWavefunction ref_wfn, Options &options) {
    auto ciwfn = std::make_shared<CIWavefunction>(ref_wfn, options);

    ciwfn->compute_energy();

    SharedWavefunction base_ciwfn = static_cast<SharedWavefunction>(ciwfn);
    return base_ciwfn;
}

}  // namespace detci
}  // namespace psi
