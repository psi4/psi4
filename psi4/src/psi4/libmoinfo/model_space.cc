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

#include "model_space.h"
#include "moinfo.h"

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cstdio>

namespace psi {

ModelSpace::ModelSpace(MOInfo* moinfo_obj_) : moinfo_obj(moinfo_obj_) {
    startup();
    build();
    classify();
    //  print();
}

ModelSpace::~ModelSpace() { cleanup(); }

void ModelSpace::startup() { wfn_sym = moinfo_obj->get_wfn_sym(); }

void ModelSpace::cleanup() {}

void ModelSpace::print() {
    outfile->Printf("\n\n  Model space:");
    outfile->Printf("\n  ------------------------------------------------------------------------------");
    for (size_t mu = 0; mu < determinants.size(); ++mu) {
        outfile->Printf("\n  %2zu %s", mu, determinants[mu].get_label().c_str());
    }
    outfile->Printf("\n\n  Closed-shell to model space mapping");
    for (size_t mu = 0; mu < closed_to_all.size(); ++mu) {
        outfile->Printf("\n  %zu -> %d", mu, closed_to_all[mu]);
    }
    outfile->Printf("\n\n  Open-shell to model space mapping");
    for (size_t mu = 0; mu < opensh_to_all.size(); ++mu) {
        outfile->Printf("\n  %zu -> %d", mu, opensh_to_all[mu]);
    }
}

}  // namespace psi
