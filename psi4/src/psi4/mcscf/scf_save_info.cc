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

#include <utility>
#include <algorithm>
#include <cstdio>

#include <cmath>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/factory.h"
#include "psi4/psifiles.h"

#include "psi4/psi4-dec.h"

#include "scf.h"

namespace psi{ namespace mcscf{

using namespace std;

void SCF::save_info()
{
    // No projection of MOs, just yet
    nmo_ = nso_;

    // figure out how many frozen orbitals per irrep
    int nfrzc = molecule_->nfrozen_core();
    intvec frz;
    for(int h = 0; h < nirreps; ++h) frz.push_back(0);
    vector<std::pair<double, int> > sorted_evals;
    for(int h = 0; h < nirreps; ++h)
      for(int i = 0; i < sopi[h]; ++i)
        sorted_evals.push_back( make_pair(epsilon->get(h,i),h) );
    sort(sorted_evals.begin(),sorted_evals.end());
    for(int i = 0; i < nfrzc; ++i)
      frz[sorted_evals[i].second]++;

    for(int h = 0; h < nirreps; ++h){
        doccpi_[h] = docc[h];
        soccpi_[h] = actv[h];
        nmopi_[h]  = nsopi_[h];
        frzcpi_[h] = frz[h];
        frzvpi_[h] = 0;
    }
    // Save the eigenvectors after rotating them
    if(options_.get_double("ROTATE_MO_ANGLE") != 0.0){
        double mo_rotate_angle = options_.get_double("ROTATE_MO_ANGLE");
        int p = options_.get_int("ROTATE_MO_P") -1;  // P, Q and IRREPS are one-based
        int q = options_.get_int("ROTATE_MO_Q") -1;
        int h = options_.get_int("ROTATE_MO_IRREP") - 1;

        outfile->Printf("\n\n  Rotating MOs %d and %d of irrep %d by %lf degrees",
                        p,q,h,mo_rotate_angle);
        double angle = mo_rotate_angle * acos(-1.0) / 180.0;
        for(int i = 0; i < sopi[h]; ++i){
            double Cp = cos(angle) * C->get(h,i,p) + sin(angle) * C->get(h,i,q);
            double Cq = cos(angle) * C->get(h,i,q) - sin(angle) * C->get(h,i,p);
            C->set(h,i,p,Cp);
            C->set(h,i,q,Cq);
        }
    }

    // Store the information in wavefunction's objects, for later access
    Ca_ = SharedMatrix(factory_->create_matrix("C"));
    Cb_ = Ca_;
    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = epsilon_a_;
    for(int h = 0; h < nirreps; ++h){
        for(int so = 0; so < nsopi_[h]; ++so){
            epsilon_a_->set(h, so, epsilon->get(h, so));
            for(int mo = 0; mo < nmopi_[h]; ++mo){
                Ca_->set(h, so, mo, C->get(h, so, mo));
            }
        }
    }

    Process::environment.globals["MCSCF TOTAL ENERGY"] = total_energy;
    Process::environment.globals["CURRENT ENERGY"] = total_energy;
    Process::environment.globals["CURRENT REFERENCE ENERGY"] = total_energy;
    energy_ = total_energy;
    cleanup();

    return;

}

}} /* End Namespaces */
