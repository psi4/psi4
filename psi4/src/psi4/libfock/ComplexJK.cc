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

#include "psi4/libfock/ComplexJK.h"

#include "psi4/pybind11.h"
// #include "psi4/libpsi4util/process.h" this also includes options as below
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/basisset.h"

namespace psi {

ComplexJK::ComplexJK(std::shared_ptr<BasisSet> primary) : BaseJK(primary) { common_init(); }

ComplexJK::~ComplexJK() {}

void ComplexJK::common_init() {
    do_csam_ = false;
    // std::shared_ptr<IntegralFactory> integral =
    //     std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
}

std::shared_ptr<ComplexJK> ComplexJK::build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                               Options& options, std::string jk_type) {

    bool do_df_scf_guess = options.get_bool("DF_SCF_GUESS");
    if (do_df_scf_guess) {
        std::string error_message = "SCREENING=DENSITY has not been implemented for ";
        error_message += (do_df_scf_guess) ? "DF_SCF_GUESS" : jk_type;
        error_message += ".";

        throw PSIEXCEPTION(error_message);
    }

    // set up ERI cutoff value
    double cutoff = options.get_str("SCREENING") == "NONE" ? 0.0 : options.get_double("INTS_TOLERANCE");

    if (jk_type == "DIRECT") {
        auto jk = std::make_shared<ComplexDirectJK>(primary, options);

        if (options["INTS_TOLERANCE"].has_changed() || options.get_str("SCREENING") == "NONE") jk->set_cutoff(cutoff);
        if (options["SCREENING"].has_changed()) jk->set_csam(options.get_str("SCREENING") == "CSAM");
        if (options["PRINT"].has_changed()) jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed()) jk->set_debug(options.get_int("DEBUG"));

        return jk;
    } else {
        throw PSIEXCEPTION("ComplexJK::build_JK: Unknown SCF Type '" + jk_type + "'");
    }
}

std::shared_ptr<ComplexJK> ComplexJK::build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                               Options& options) {
    return build_JK(primary, auxiliary, options, options.get_str("SCF_TYPE"));
}

std::shared_ptr<ComplexJK> ComplexJK::build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                               Options& options, bool /*do_wK*/, size_t /*doubles*/) {
    // This constructor still exists to reflect the same interface as JK::build_JK
    return build_JK(primary, auxiliary, options, options.get_str("SCF_TYPE"));
}

void ComplexJK::allocate_JK() {
    if (J_.size() != D_.size()) {
        J_.clear();
        K_.clear();

        for (size_t i = 0; i < D_.size() && do_J_; i++) {
            std::stringstream s;
            s << "J " << i << " (SO)";
            // J_.push_back(std::make_shared<ComplexMatrix>(s.str(), D_[i]->rowspi(), D_[i]->rowspi(),
            //                                       D_[i]->symmetry()));
        }

    }
}

void ComplexJK::compute_D() { throw pybind11::attribute_error("compute_D not implemented!"); }

size_t ComplexJK::memory_overhead() const {
    size_t mem = 0L;
    int JKwKD_factor = 1;
    if (do_J_) JKwKD_factor++;
    if (do_K_) JKwKD_factor++;

    int C_factor = 1;
    if (!lr_symmetric_) C_factor++;

    // if (C1() && C_left_.size() && C_left_[0]->num_blocks() != 1) {
    //     int nbf = primary_->nbf();
    //     for (size_t N = 0; N < C_left_.size(); N++) {
    //         int nocc = 0;
    //         for (int h = 0; h < C_left_[N]->num_blocks(); h++) {
    //             nocc += C_left_[N]->row_dims()[h];
    //         }
    //         mem += C_factor * (size_t)nocc * nbf + JKwKD_factor * (size_t)nbf * nbf;
    //     }
    // }

    return mem;
}

void ComplexJK::zero() {
    if (do_J_) {
        for(auto J : J_) J->zero();
    }
    if (do_K_) {
        for(auto K : K_) K->zero();
    }
}

void ComplexJK::compute() { throw pybind11::attribute_error("ComplexJK::compute not implemented!"); }

}  // namespace psi
