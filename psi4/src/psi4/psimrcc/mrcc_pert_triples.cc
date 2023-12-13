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

/**
 *  @file ccmrcc_pert_triples.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Computes the (T) correction
 */

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/process.h"

#include "psi4/pragma.h"
#include <memory>
#include "psi4/psi4-dec.h"

#include "mrcc.h"
#include "mrccsd_t.h"

namespace psi {
namespace psimrcc {

void CCMRCC::compute_perturbative_triples() {
    Timer timer;

    h_eff.set_eigenvalue(current_energy);
    h_eff.set_matrix(Heff, wfn_->moinfo()->get_nrefs());
    h_eff.set_right_eigenvector(right_eigenvector);
    h_eff.set_left_eigenvector(left_eigenvector);
    h_eff.set_zeroth_order_eigenvector(zeroth_order_eigenvector);

    MRCCSD_T mrccsd_t(options_, &h_eff);

    if (options_.get_bool("DIAGONALIZE_HEFF")) {
        outfile->Printf("\n\n  Diagonalizing Heff");
        current_energy = h_eff.diagonalize();
    } else {
        outfile->Printf("\n\n  Computing the expectation value of Heff");
        current_energy = h_eff.expectation_value();
    }
    wfn_->set_scalar_variable("CURRENT ENERGY", current_energy);
    wfn_->set_scalar_variable("MRCC TOTAL ENERGY", current_energy);

    outfile->Printf("\n\n%6c* Mk-MRCCSD(T) total energy   =    %20.12f", ' ', current_energy);
    outfile->Printf("\n\n  Timing for triples:             %20.6f s", timer.get());
}

}  // namespace psimrcc
}  // namespace psi
