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
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/** Standard library includes */
#include <fstream>
#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/sieve.h"
#include "dfocc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::dfgrad()
{

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
    }
    else {
        tpdm_tilde_cc();
        back_trans_cc();
    }

//===========================================================================================
//============================ Gradient =====================================================
//===========================================================================================
    outfile->Printf("\tComputing analytic gradients...\n");


    gradient_terms.push_back("Nuclear");
    gradient_terms.push_back("Kinetic");
    gradient_terms.push_back("Potential");
    gradient_terms.push_back("Overlap");
    gradient_terms.push_back("3-Index:RefSep");
    gradient_terms.push_back("3-Index:Corr");
    gradient_terms.push_back("Metric:RefSep");
    gradient_terms.push_back("Metric:Corr");
    gradient_terms.push_back("Total");

    // OEI GRAD
    oei_grad();

    // TEI GRAD
    tei_grad_ref();
    tei_grad_corr();

//===========================================================================================
//========================= Total Gradient ==================================================
//===========================================================================================
    // => Total Gradient <= //
    SharedMatrix total = SharedMatrix(gradients["Nuclear"]->clone());
    total->zero();

    for (int i = 0; i < gradient_terms.size(); i++) {
        if (gradients.count(gradient_terms[i])) {
            total->add(gradients[gradient_terms[i]]);
        }
    }

    gradients["Total"] = total;
    gradients["Total"]->set_name("Total Gradient");

    // OEI grad
    gradients["One-Electron"] = SharedMatrix(gradients["Nuclear"]->clone());
    gradients["One-Electron"]->set_name("One-Electron Gradient");
    gradients["One-Electron"]->zero();
    gradients["One-Electron"]->add(gradients["Kinetic"]);
    gradients["One-Electron"]->add(gradients["Potential"]);
    gradients["One-Electron"]->print_atom_vector();

    // TEI grad
    gradients["Two-Electron"] = SharedMatrix(gradients["Nuclear"]->clone());
    gradients["Two-Electron"]->set_name("Two-Electron Gradient");
    gradients["Two-Electron"]->zero();
    gradients["Two-Electron"]->add(gradients["3-Index:RefSep"]);
    gradients["Two-Electron"]->add(gradients["3-Index:Corr"]);
    gradients["Two-Electron"]->add(gradients["Metric:RefSep"]);
    gradients["Two-Electron"]->add(gradients["Metric:Corr"]);
    gradients["Two-Electron"]->print_atom_vector();//UB


    // => Final Printing <= //
    if (print_ > 1) {
        for (int i = 0; i < gradient_terms.size(); i++) {
            if (gradients.count(gradient_terms[i])) {
                gradients[gradient_terms[i]]->print_atom_vector();
            }
        }
    } else {
        gradients["Total"]->print_atom_vector();
    }

    gradient_ = total;

//outfile->Printf("\tdfgrad is done. \n");
}// end dfgrad

}} // End Namespaces


