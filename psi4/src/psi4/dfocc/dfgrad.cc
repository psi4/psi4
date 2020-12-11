/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

/** Standard library includes */
#include <fstream>
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "dfocc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::dfgrad() {
    //===========================================================================================
    //============================ Preliminaries ================================================
    //===========================================================================================
    tstop();
    tstart();
    title_grad();
    outfile->Printf("\tAnalytic gradients computation is starting...\n");

    if (wfn_type_ == "DF-OMP2") {
        tpdm_tilde();
        back_trans();
    } else {
        tpdm_tilde_cc();
        back_trans_cc();
    }

    //===========================================================================================
    //============================ Gradient =====================================================
    //===========================================================================================
    outfile->Printf("\tComputing analytic gradients...\n");

    std::map<std::string, SharedMatrix> gradients;

    // OEI GRAD
    oei_grad(gradients);

    // TEI GRAD
    tei_grad("JK", gradients);
    tei_grad("RI", gradients);

    //===========================================================================================
    //========================= Total Gradient ==================================================
    //===========================================================================================
    // => Total Gradient <= //
    auto total = SharedMatrix(gradients["Nuclear"]->clone());
    total->zero();

    for (auto& kv: gradients) {
        total->add(kv.second);
    }

    // OEI grad
    gradients["One-Electron"] = gradients["Core"]->clone();
    gradients["One-Electron"]->set_name("One-Electron Gradient");
    gradients["One-Electron"]->print_atom_vector();

    // TEI grad
    gradients["Two-Electron"] = SharedMatrix(gradients["Nuclear"]->clone());
    gradients["Two-Electron"]->set_name("Two-Electron Gradient");
    gradients["Two-Electron"]->zero();
    gradients["Two-Electron"]->add(gradients["3-Index:RefSep"]);
    gradients["Two-Electron"]->add(gradients["3-Index:Corr"]);
    gradients["Two-Electron"]->add(gradients["Metric:RefSep"]);
    gradients["Two-Electron"]->add(gradients["Metric:Corr"]);
    gradients["Two-Electron"]->print_atom_vector();  // UB

    gradients["Total"] = total;
    gradients["Total"]->set_name("Total Gradient");

    // => Final Printing <= //
    if (print_ > 1) {
        for (auto& kv: gradients) {
            kv.second->print_atom_vector();
        }
    } else {
        gradients["Total"]->print_atom_vector();
    }

    set_gradient(total);

    // outfile->Printf("\tdfgrad is done. \n");
}  // end dfgrad

}  // namespace dfoccwave
}  // namespace psi
