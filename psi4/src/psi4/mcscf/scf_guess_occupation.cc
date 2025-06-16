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

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <cstdio>

#include "psi4/libmoinfo/libmoinfo.h"

#include "scf.h"

extern FILE* outfile;

namespace psi {
namespace mcscf {

void SCF::guess_occupation() {
    if (moinfo_scf->get_guess_occupation()) {
        // Assumes the eigenvalues of some Fock operator
        // are in the SBlockVector epsilon
        std::vector<std::pair<double, int> > evals;

        for (int h = 0; h < nirreps; ++h)
            for (int i = 0; i < sopi[h]; ++i) evals.push_back(std::make_pair(epsilon->get(h, i), h));

        // Sort the eigenvalues by energy
        sort(evals.begin(), evals.end());

        int ndocc = std::min(moinfo_scf->get_nael(), moinfo_scf->get_nbel()) - (reference == tcscf ? 1 : 0);
        int nactv = std::abs(moinfo_scf->get_nael() - moinfo_scf->get_nbel()) + (reference == tcscf ? 2 : 0);

        std::vector<int> new_docc;
        std::vector<int> new_actv;
        for (int h = 0; h < nirreps; ++h) {
            new_docc.push_back(0);
            new_actv.push_back(0);
        }
        for (int i = 0; i < ndocc; ++i) {
            new_docc[evals[i].second]++;
        }
        for (int i = ndocc; i < ndocc + nactv; ++i) {
            new_actv[evals[i].second]++;
        }

        if ((new_docc != docc) || (new_actv != actv)) {
            outfile->Printf("\n\n  Occupation changed");
            outfile->Printf("\n  docc = (");
            for (int h = 0; h < nirreps; ++h) outfile->Printf(" %d", new_docc[h]);
            outfile->Printf(" )");
            outfile->Printf("\n  actv = (");
            for (int h = 0; h < nirreps; ++h) outfile->Printf(" %d", new_actv[h]);
            outfile->Printf(" )\n");
        }
        docc = new_docc;
        actv = new_actv;
    }
}

}  // namespace mcscf
}  // namespace psi
