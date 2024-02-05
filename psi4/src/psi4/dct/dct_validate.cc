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

#include "dct.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace dct {

void DCTSolver::validate_energy() {
    if (options_.get_str("ALGORITHM") == "TWOSTEP" && orbital_optimized_)
        throw PSIEXCEPTION("Two-step algorithm cannot be run for orbital-optimized DCT methods");

    // RHF-specific warnings
    if (same_a_b_orbs_) {
        const std::string msg = "RHF-reference DCT";
        if (options_.get_str("THREE_PARTICLE") == "PERTURBATIVE")
            throw FeatureNotImplemented(msg, "Three-particle energy correction", __FILE__, __LINE__);
        if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13")
            throw FeatureNotImplemented(msg, "DCT_FUNCTIONAL = ODC-13", __FILE__, __LINE__);
        if (options_.get_str("ALGORITHM") == "QC" || options_.get_str("ALGORITHM") == "TWOSTEP")
            throw FeatureNotImplemented(msg, "ALGORITHM = QC/TWOSTEP", __FILE__, __LINE__);
        if (options_.get_str("DCT_GUESS") == "CC" || options_.get_str("DCT_GUESS") == "BCC") {
            throw FeatureNotImplemented(msg, "CC/BCC amplitude guess", __FILE__, __LINE__);
        }
    }

    if (options_.get_str("DCT_TYPE") == "DF") {
        const std::string msg = "Density-fitted DCT";
        if (options_.get_str("THREE_PARTICLE") == "PERTURBATIVE")
            throw FeatureNotImplemented(msg, "Three-particle energy correction", __FILE__, __LINE__);
        if (options_.get_str("DCT_FUNCTIONAL") == "ODC-13")
            throw FeatureNotImplemented(msg, "DCT_FUNCTIONAL = ODC-13", __FILE__, __LINE__);
        if (options_.get_str("ALGORITHM") == "QC")
            throw FeatureNotImplemented(msg, "ALGORITHM = QC", __FILE__, __LINE__);
    }

    if (options_.get_str("ALGORITHM") == "QC") {
        if (options_.get_str("AO_BASIS") == "DISK" && options_.get_str("QC_TYPE") == "SIMULTANEOUS")
            throw FeatureNotImplemented("Simultaneous QC", "AO_BASIS = DISK", __FILE__, __LINE__);
        if (options_.get_str("DCT_FUNCTIONAL") != "DC-06")
            outfile->Printf(
                "\n\n\t**** Warning: Using DC-06 hessian, as others not implemented. Quadratic convergence is not "
                "guaranteed. ****\n");
    }

    if ((options_.get_str("SCF_TYPE").find("DF") != std::string::npos) || options_.get_str("SCF_TYPE") == "CD" ||
        options_.get_str("SCF_TYPE") == "DIRECT") {
        if (!options_["DCT_TYPE"].has_changed())
            options_.set_global_str("DCT_TYPE", "DF");
        else if (options_.get_str("DCT_TYPE") == "CONV")
            throw PSIEXCEPTION("Please set SCF_TYPE to PK or OUT_OF_CORE in order to use DCT_TYPE=CONV.");
    }

    if (options_.get_str("DCT_TYPE") == "DF") {
        if (!options_["AO_BASIS"].has_changed())
            options_.set_str("DCT", "AO_BASIS", "NONE");
        else if (options_.get_str("AO_BASIS") == "DISK") {
            outfile->Printf(
                "\n\n\t**** Warning: AO_BASIS=DISK not implemented in density-fitted DCT. Switching to AO_BASIS=NONE "
                "****\n");
            options_.set_str("DCT", "AO_BASIS", "NONE");
        }
    }

    if (options_.get_str("DCT_FUNCTIONAL") == "CEPA0")
        throw PSIEXCEPTION(
            "CEPA0 was removed from the DCT module in 1.4. Please use the lccd method in OCC, DFOCC, or FNOCC.");
}

void DCTSolver::validate_opdm() {
    if (!options_.get_bool("OPDM")) return;
    const std::string msg = "Density matrices, properties, and analytic gradients";
    if (options_.get_str("DCT_FUNCTIONAL") == "DC-12")
        throw FeatureNotImplemented("DC-12 functional", msg, __FILE__, __LINE__);
    if (options_.get_str("THREE_PARTICLE") == "PERTURBATIVE")
        throw FeatureNotImplemented("Three-particle energy correction", msg, __FILE__, __LINE__);
    if (options_.get_str("DCT_FUNCTIONAL") == "DC-06" && same_a_b_orbs_)
        throw FeatureNotImplemented("RHF DC-06", msg, __FILE__, __LINE__);
}

// At present, no special validation tasks needed.
void DCTSolver::validate_gradient() {}

}  // namespace dct
}  // namespace psi
