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
    // Not sure if I need this yet
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
        throw PSIEXCEPTION("ComplexJK::build_JK: Unknown complex SCF Type '" + jk_type + "'");
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

        for (size_t N = 0; N < D_.size() && do_J_; N++) {
            std::stringstream s;
            s << "J " << N << " (SO)";
            J_.push_back(std::make_shared<ComplexMatrix>(s.str(), D_[N]->block_dims()));
        }

        for (size_t N = 0; N < D_.size() && do_K_; N++) {
            std::stringstream s;
            s << "K " << N << " (SO)";
            K_.push_back(std::make_shared<ComplexMatrix>(s.str(), D_[N]->block_dims()));
        }
    }
}

void ComplexJK::compute_D() {
    /// Make sure the memory is there
    bool same = true;
    if (C_left_.size() != D_.size()) {
        same = false;
    }

    if (!same) {
        D_.clear();
        for (size_t N = 0; N < C_left_.size(); N++) {
            std::stringstream s;
            s << "D " << N << " (SO)";
            // eventually maybe switch to TiledTensor so C_left_ and C_right_ can
            // have different number of rows.
            D_.push_back(std::make_shared<ComplexMatrix>(s.str(), C_left_[N]->block_dims()));
        }
    }

    // Form the density, differs from dou
    for (size_t N = 0; N < D_.size(); ++N) {
        // int symm = D_[N]->symmetry();
        D_[N]->zero();
        auto temp1_ = ComplexMatrix("square C", D_[N]->block_dims());
        auto temp2_ = ComplexMatrix("square C.H", D_[N]->block_dims());

        for (size_t h = 0; h < D_[N]->num_blocks(); h++) {
            // C_DGEMM('N', 'T', nsol, nsor, nocc, 1.0, Clp[0], nocc, Crp[0], nocc, 0.0, Dp[0], nsor);

            temp1_[h].zero();
            temp2_[h].zero();

            /// Yeah, this one isn't ready yet.
            throw pybind11::attribute_error("DirectJK::compute_D() not implemented!");

            // Fills temp1_ and temp2_ with the occupied (2*nsopi_ x nelecpi_) and conjugate
            // occupied matrices (e.g. Cocc)
            // Note that temp1_ and temp2_ are both (2*nsopi_ x 2*nsopi_), but the 'remainder'
            // are all->zeros and contribute nothing
            // This is just a minor inefficiency, since these matrices are larger than they should be
            // NOTE: this is done like this because BlockTensors cannot handle rectangular matrices

            // TODO: get nelecpi_[h] somehow
            auto nelecpi_ = C_left_[N]->block_dim(h);

            for (int j = 0; j < C_left_[N]->block_dim(h); j++) {
                for (int k = 0; k < nelecpi_; k++) {
                    const std::complex<double>& C_jk = C_left_[N]->block(h)(j, k);
                    temp1_[h](j, k) = C_jk;
                    temp2_[h](j, k) = std::conj(C_jk);
                }
            }
        }

        // D_ = einsums("ui,vi->uv", temp1_, temp2_)
        einsums::tensor_algebra::einsum(einsums::Indices{einsums::index::u, einsums::index::v}, D_[N].get(),  // D_uv
                                        einsums::Indices{einsums::index::u, einsums::index::i}, temp1_,   // Cocc_ui
                                        einsums::Indices{einsums::index::v, einsums::index::i}, temp2_    // Cocc_vi.conj().T
        );

    }
}

size_t ComplexJK::memory_overhead() const {
    // Only non-zero for USO / AO quantities
    return 0L;
}

void ComplexJK::zero() {
    if (do_J_) {
        for(auto J : J_) J->zero();
    }
    if (do_K_) {
        for(auto K : K_) K->zero();
    }
}

void ComplexJK::compute() {
    // Is this density symmetric?
    if (C_left_.size() && !C_right_.size()) {
        lr_symmetric_ = true;
        C_right_ = C_left_;
    } else {
        throw PSIEXCEPTION("C_left must equal C_right.");
    }

    // input_symmetry_cast_map_ not needed for C1

    // Construct the densities
    timer_on("JK: D");
    compute_D();
    timer_off("JK: D");

    allocate_JK();

    timer_on("JK: JK");
    compute_JK();
    timer_off("JK: JK");

    if (debug_ > 6) {
        outfile->Printf("\n   WARNING: ComplexMatrix->print() not implemented. "
                        "No internal JK quantities will be logged.\n");
        // for (size_t N = 0; N < C_left_.size(); N++) {
        //     C_left_[N]->print("outfile");
        //     C_right_[N]->print("outfile");
        //     D_[N]->print("outfile");
        //     J_[N]->print("outfile");
        //     K_[N]->print("outfile");
        // }
    }

    // Assuming lr_symmetric_
    C_right_.clear();
}

}  // namespace psi
