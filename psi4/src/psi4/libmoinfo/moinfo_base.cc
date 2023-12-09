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

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "psi4/psi4-dec.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "moinfo_base.h"

namespace psi {

/// @brief Constructor for MOInfoBase class
/// @param ref_wfn_
/// @param options_
/// @param silent_
MOInfoBase::MOInfoBase(Wavefunction& ref_wfn_, Options& options_, bool silent_)
    : options(options_),
      silent(silent_),
      ref_wfn(ref_wfn_),
      nirreps(ref_wfn_.nirrep()),
      charge(ref_wfn_.molecule()->molecular_charge()) {
    nso = 0;
    nmo = 0;
    ndocc = 0;
    nactv = 0;
    nael = 0;
    nbel = 0;
    nactive_ael = 0;
    nactive_bel = 0;
    wfn_sym = 0;
    guess_occupation = true;

    multiplicity = ref_wfn.molecule()->multiplicity();
}

void MOInfoBase::read_data() {
    nso = ref_wfn.nso();
    // Read sopi and save as a STL vector
    if (nirreps != ref_wfn.nsopi().n())
        throw PSIEXCEPTION(
            "\n\n  MOInfoBase::read_data(): Suspicious condition! The number of irreps in the reference wavefunction "
            "is not equal to the size of the number of SOs per irrep array.\n\n");
    sopi = ref_wfn.nsopi().blocks();
    irr_labs = ref_wfn.molecule()->irrep_labels();
    nuclear_energy = ref_wfn.molecule()->nuclear_repulsion_energy(ref_wfn.get_dipole_field_strength());
}

void MOInfoBase::compute_number_of_electrons() {
    int nel = 0;
    int natom = ref_wfn.molecule()->natom();

    for (int i = 0; i < natom; i++) {
        nel += static_cast<int>(ref_wfn.molecule()->Z(i));
    }
    nel -= charge;

    // Check if the multiplicity makes sense
    if (((nel + 1 - multiplicity) % 2) != 0) throw PSIEXCEPTION("\n\n  MOInfoBase: Wrong multiplicity.\n\n");
    nael = (nel + multiplicity - 1) / 2;
    nbel = nel - nael;
}

void MOInfoBase::read_mo_space(int nirreps_ref, int& n, intvec& mo, std::string labels) {
    bool read = false;

    std::vector<std::string> label_vec = split(labels);
    for (size_t k = 0; k < label_vec.size(); ++k) {
        // Does the array exist in the input?
        std::string& label = label_vec[k];
        if (!options[label].has_changed()) continue;  // The user didn't specify this, it's just the default
        int size = options[label].size();
        // Defaults is to set all to zero
        mo.assign(nirreps_ref, 0);
        n = 0;
        if (read) {
            outfile->Printf("\n\n  libmoinfo has found a redundancy in the input keywords %s , please fix it!",
                            labels.c_str());

            exit(1);
        } else {
            read = true;
        }
        if (size == nirreps_ref) {
            for (int i = 0; i < size; i++) {
                mo[i] = options[label][i].to_integer();
                n += mo[i];
            }
        } else {
            outfile->Printf(
                "\n\n  The size of the %s array (%d) does not match the number of irreps (%d), please fix the input "
                "file",
                label_vec[k].c_str(), size, nirreps_ref);

            exit(1);
        }
    }
}

void MOInfoBase::print_mo_space(int n, intvec& mo, std::string labels) {
    outfile->Printf("\n  %s", labels.c_str());

    for (int i = nirreps; i < 8; i++) outfile->Printf("     ");
    for (int i = 0; i < nirreps; i++) outfile->Printf(" %3d ", mo[i]);
    outfile->Printf("  %3d", n);
}

intvec MOInfoBase::convert_int_array_to_vector(int n, const int* array) {
    // Read an integer array and save as a STL vector
    return intvec(array, array + n);
}

}  // namespace psi
