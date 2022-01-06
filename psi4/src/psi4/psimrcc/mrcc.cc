/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "algebra_interface.h"
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "updater.h"

namespace psi {
namespace psimrcc {

CCMRCC::CCMRCC(std::shared_ptr<PSIMRCCWfn> wfn, Options &options)
    : CCManyBody(wfn, options), options_(options), h_eff(Hamiltonian(wfn)) {
    triples_type = ccsd;
    triples_coupling_type = cubic;
    ap_correction = false;  // Set tu true when computing the a posteriori correction
    current_energy = 0.0;
    old_energy = 10.0;

    // Parse the CORR_WFN parameter
    std::vector<std::string> theory_levels = split("PT2 CCSD CCSD_T CCSDT-1A CCSDT-1B CCSDT-2 CCSDT-3 CCSDT");
    for (size_t i = 0; i < theory_levels.size(); ++i) {
        if (options.get_str("CORR_WFN") == theory_levels[i]) triples_type = TriplesType(i);
    }

    // Parse the COUPLING parameter
    std::vector<std::string> coupling_levels = split("NONE LINEAR QUADRATIC CUBIC");
    for (size_t i = 0; i < coupling_levels.size(); ++i) {
        if (options.get_str("COUPLING") == coupling_levels[i]) {
            triples_coupling_type = TriplesCouplingType(i);
        }
    }

    // Add the matrices that will store the intermediates
    add_matrices();

    // Generate the Fock matrices, Integrals and Denominators
    generate_integrals();
    generate_denominators();

    if (triples_type > ccsd) generate_triples_denominators();

    compute_reference_energy();

    // Initialize the appropriate updater
    if (options.get_str("CORR_ANSATZ") == "MK")
        updater_ = std::make_shared<MkUpdater>(wfn, options);
    else if (options.get_str("CORR_ANSATZ") == "BW")
        updater_ = std::make_shared<BWUpdater>(wfn, options);
    else
        throw PSIEXCEPTION("I don't know what updater goes with CORR_ANSATZ " + options.get_str("CORR_ANSATZ"));
}

CCMRCC::~CCMRCC() {}

}  // namespace psimrcc
}  // namespace psi
