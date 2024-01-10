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

#include "dct.h"
#include "psi4/psifiles.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.h"

namespace psi {
namespace dct {

SharedMatrix DCTSolver::compute_gradient() {
    // Print out the header
    outfile->Printf("\n\n\t***********************************************************************************\n");
    outfile->Printf("\t*                            DCT Analytic Gradients Code                          *\n");
    outfile->Printf("\t*                by Alexander Sokolov, Andy Simmonett, and Xiao Wang              *\n");
    outfile->Printf("\t***********************************************************************************\n\n");

    validate_energy();
    validate_opdm();
    validate_gradient();
    if (orbital_optimized_) oo_gradient_init();

    // If the system is closed-shell, then ...
    if (options_.get_str("REFERENCE") == "RHF") {
        compute_gradient_RHF();
    }
    // If the system is open-shell, then ...
    else {
        compute_gradient_UHF();
    }

    return std::make_shared<Matrix>("nullptr", 0, 0);
}
}  // namespace dct
}  // namespace psi
