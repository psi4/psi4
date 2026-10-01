/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "mp2.h"
#include "streamed.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"

namespace psi {
namespace f12 {

SharedWavefunction f12(SharedWavefunction ref_wfn, Options& options) {
    std::shared_ptr<Wavefunction> f12;
    if (options.get_str("F12_SUBTYPE") == "STREAMED") {
        if (options.get_str("MP2_TYPE") != "DF") {
            throw PSIEXCEPTION("F12_SUBTYPE=STREAMED requires MP2_TYPE=DF");
        }
        if (options.get_bool("F12_READ_INTS")) {
            throw PSIEXCEPTION("F12_SUBTYPE=STREAMED does not read saved F12 integrals");
        }
        if (options.get_int("F12_AUX_BLOCK_SIZE") < 1) {
            throw PSIEXCEPTION("F12_AUX_BLOCK_SIZE must be positive");
        }
        if (ref_wfn->molecule()->schoenflies_symbol() != "c1") {
            throw PSIEXCEPTION("F12_SUBTYPE=STREAMED requires C1 symmetry");
        }
        if (ref_wfn->basisset()->has_ECP() || ref_wfn->get_basisset("CABS")->has_ECP() ||
            ref_wfn->get_basisset("DF_BASIS_MP2")->has_ECP()) {
            throw PSIEXCEPTION("F12_SUBTYPE=STREAMED currently supports all-electron basis sets only");
        }
        if (ref_wfn->frzvpi()[0] != 0 || ref_wfn->doccpi()[0] <= ref_wfn->frzcpi()[0] ||
            ref_wfn->nmo() <= ref_wfn->doccpi()[0]) {
            throw PSIEXCEPTION("F12_SUBTYPE=STREAMED requires active occupied and virtual orbitals, with no frozen virtuals");
        }
        f12 = std::make_shared<StreamedMP2F12>(ref_wfn, options);
    } else if (options.get_str("F12_SUBTYPE").find("DISK") != std::string::npos) {
        f12 = std::make_shared<DiskMP2F12>(ref_wfn, options);
    } else {
        f12 = std::make_shared<MP2F12>(ref_wfn, options);
    }

    return f12;
}

}  // namespace f12
}  // namespace psi
