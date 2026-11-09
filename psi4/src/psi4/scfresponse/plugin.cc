/*
 * @BEGIN LICENSE
 *
 * scf_response by Psi4 Developer, a plugin to:
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

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.hpp"

#include "uks_response.h"
#include "rks_response.h"

namespace psi{ namespace scf_response {

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    if (name == "SCF_RESPONSE"|| options.read_globals()) {
        /*- response properties -*/
        options.add_str("PROPERTY", "POLARIZABILITY", "POLARIZABILITY HYPERPOLARIZABILITY SHG OR POCKELS ALL");
    }

    return true;
}

extern "C" PSI_API
SharedWavefunction scf_response(SharedWavefunction ref_wfn, Options& options)
{
    // response properties
    if ( options["PROPERTY"].has_changed() ) {
        if ( options.get_str("REFERENCE") == "UKS" || options.get_str("REFERENCE") == "UHF" ) {
            std::shared_ptr<UKSResponseSolver> uks (new UKSResponseSolver((std::shared_ptr<Wavefunction>)ref_wfn, options));
            uks->compute_properties();
        }else if ( options.get_str("REFERENCE") == "RHF" || options.get_str("REFERENCE") == "RKS" ) {
            std::shared_ptr<RKSResponseSolver> rks (new RKSResponseSolver((std::shared_ptr<Wavefunction>)ref_wfn, options));
            rks->compute_properties();
        }
    }else {
        throw PsiException("hmm probably you meant to set PROPERTY",__FILE__,__LINE__);
    }
    
    return ref_wfn;
}

}} // End namespaces

