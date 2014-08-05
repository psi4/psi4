/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/** Standard library includes */
#include <fstream>
#include <psifiles.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <libmints/sieve.h>
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
    title_grad();
    fprintf(outfile,"\tAnalytic gradients computation is starting...\n");
    fflush(outfile);
    tpdm_tilde();
    back_trans();

//===========================================================================================
//============================ Gradient =====================================================
//===========================================================================================
    fprintf(outfile,"\tComputing analytic gradients...\n");  
    fflush(outfile);

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

    Process::environment.set_gradient(total);
    Process::environment.wavefunction()->set_gradient(total);


//fprintf(outfile,"\tdfgrad is done. \n"); fflush(outfile);
}// end dfgrad 

}} // End Namespaces




