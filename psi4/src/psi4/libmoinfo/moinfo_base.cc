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

MOInfoBase::MOInfoBase(Wavefunction& ref_wfn_, Options& options_)
    : options(options_), ref_wfn(ref_wfn_), charge(ref_wfn_.molecule()->molecular_charge()) {
    nso = 0;
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
    nirreps = ref_wfn.nirrep();
    nso = ref_wfn.nso();
    // Read sopi and save as a STL vector
    if (ref_wfn.nirrep() != ref_wfn.nsopi().n()) {
        const std::string msg =
            "MOInfoBase::read_data(): Suspicious condition! The number of irreps in the reference wavefunction is not "
            "equal to the size of the number of SOs per irrep array. Invalid ref_wfn?\n";
        outfile->Printf(msg.c_str());
        throw PSIEXCEPTION(msg);
    }
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

void MOInfoBase::read_mo_space(const int nirreps_ref, int& n, intvec& mo, const std::string& labels) {
    bool read = false;

    const std::vector<std::string> label_vec = split(labels);
    for (size_t k = 0; k < label_vec.size(); ++k) {
        // Does the array exist in the input?
        const std::string& label = label_vec[k];
        if (!options[label].has_changed()) continue;  // The user didn't specify this, it's just the default
        int size = options[label].size();
        // Defaults is to set all to zero
        mo.assign(nirreps_ref, 0);
        n = 0;
        if (read) {
            throw std::runtime_error("libmoinfo has found a redundancy in the input keywords " + labels +
                                     ", please fix it!\n");
        } else {
            read = true;
        }
        if (size == nirreps_ref) {
            for (int i = 0; i < size; i++) {
                mo[i] = options[label][i].to_integer();
                n += mo[i];
            }
        } else {
            std::ostringstream oss;
            oss << "The size of the " << label_vec[k] << " array (" << size << ") does not match the number of irreps ("
                << nirreps_ref << "), please fix the input\n";
            throw std::runtime_error(oss.str());
        }
    }
}

void MOInfoBase::print_mo_space(int n, const intvec& mo, const std::string& labels) {
    outfile->Printf("\n  %s", labels.c_str());

    for (int i = nirreps; i < 8; i++) outfile->Printf("     ");
    for (int i = 0; i < nirreps; i++) outfile->Printf(" %3d ", mo[i]);
    outfile->Printf("  %3d", n);
}

}  // namespace psi
