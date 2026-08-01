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

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libqt/qt.h"

#include <Einsums/BLAS.hpp>
#include <Einsums/Print.hpp>

#include <complex>
#include <sstream>

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
        error_message += ". Set DF_SCF_GUESS=False";

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
    } else if (jk_type == "PK" || jk_type == "DF") {
        throw PSIEXCEPTION("SCF type '" + jk_type + "' not implemented for complex wavefunctions. "
                           "Use SCF_TYPE=DIRECT instead.");
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
            J_.push_back(std::make_shared<ComplexMatrix>(s.str(), D_[N]->tile_size(0)));
        }

        for (size_t N = 0; N < D_.size() && do_K_; N++) {
            std::stringstream s;
            s << "K " << N << " (SO)";
            K_.push_back(std::make_shared<ComplexMatrix>(s.str(), D_[N]->tile_size(0)));
        }
    }
}

void ComplexJK::compute_D() {
    // Ensure D_ matches C_left_/C_right_ (AO_left x AO_right per irrep; square when bases match).
    if (C_left_.size() != D_.size()) {
        D_.clear();
        for (size_t N = 0; N < C_left_.size(); N++) {
            std::stringstream s;
            s << "D " << N << " (SO)";
            auto const& row_sizes = C_left_[N]->tile_size(0);
            auto const& col_sizes = (N < C_right_.size()) ? C_right_[N]->tile_size(0) : row_sizes;
            D_.push_back(std::make_shared<ComplexMatrix>(s.str(), row_sizes, col_sizes));
        }
    }

    for (size_t N = 0; N < D_.size(); ++N) {
        D_[N]->zero();

        auto const& Cl = *C_left_[N];
        // compute() aliases C_right_ = C_left_ when only C_left is filled; allow the same here.
        auto const& Cr = (N < C_right_.size()) ? *C_right_[N] : Cl;

        if (Cl.grid_size(0) != Cr.grid_size(0) || Cl.grid_size(1) != Cr.grid_size(1)) {
            throw PSIEXCEPTION("ComplexJK::compute_D: C_left/C_right tile grids must match.");
        }

        for (int h = 0; h < static_cast<int>(Cl.grid_size(0)); ++h) {
            int nsol = Cl.tile_size(0)[h];
            int nsor = Cr.tile_size(0)[h];
            int nocc = Cl.tile_size(1)[h];
            if (nocc != Cr.tile_size(1)[h]) {
                throw PSIEXCEPTION("ComplexJK::compute_D: C_left/C_right occupied dimensions must match.");
            }
            if (!nsol || !nsor || !nocc) continue;
            if (!Cl.has_tile(h, h) || !Cr.has_tile(h, h)) continue;

            auto const& Cl_h = Cl.tile(h, h);
            auto const& Cr_h = Cr.tile(h, h);
            auto& D_h = D_[N]->tile(h, h);  // allocates (nsol x nsor)

            // D_h = Cl_h * Cr_h^H
            einsums::blas::gemm('n', 'c', nsol, nsor, nocc, std::complex<double>{1.0}, Cl_h.data(),
                                static_cast<int>(Cl_h.stride(0)), Cr_h.data(), static_cast<int>(Cr_h.stride(0)),
                                std::complex<double>{0.0}, D_h.data(), static_cast<int>(D_h.stride(0)));
        }
    }
}

size_t ComplexJK::memory_overhead() const {
    // TODO: Check this?
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
        lr_symmetric_ = false;
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
        for (size_t N = 0; N < C_left_.size(); N++) {
            einsums::fprintln(*outfile->stream(), *C_left_[N]);
            einsums::fprintln(*outfile->stream(), *C_right_[N]);
            einsums::fprintln(*outfile->stream(), *D_[N]);
            if (N < J_.size()) einsums::fprintln(*outfile->stream(), *J_[N]);
            if (N < K_.size()) einsums::fprintln(*outfile->stream(), *K_[N]);
        }
    }

    if (lr_symmetric_) {
        C_right_.clear();
    }
}

}  // namespace psi
