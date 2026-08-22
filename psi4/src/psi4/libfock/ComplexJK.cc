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
#include "psi4/libmints/integral.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/matrix.h"
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

    // Set up AO2USO matrix for symmetry handling
    std::shared_ptr<IntegralFactory> integral =
        std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    auto pet = std::make_shared<PetiteList>(primary_, integral);
    AO2USO_ = std::make_shared<Matrix>(pet->aotoso());
    c1_ = (AO2USO_->nirrep() > 1);
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

        for (int N = 0; N < D_.size() && do_J_; N++) {
            std::stringstream s;
            s << "J " << N << " (SO)";
            const auto& Dt = *D_[N];
            J_.push_back(std::make_shared<ComplexMatrix>(s.str(), Dt.rowspi(), Dt.colspi()));
        }

        for (int N = 0; N < D_.size() && do_K_; N++) {
            std::stringstream s;
            s << "K " << N << " (SO)";
            const auto& Dt = *D_[N];
            K_.push_back(std::make_shared<ComplexMatrix>(s.str(), Dt.rowspi(), Dt.colspi()));
        }
    }
}

void ComplexJK::compute_D() {
    // Ensure D_ matches C_left_/C_right_ (AO_left x AO_right per irrep; square when bases match).
    if (C_left_.size() != D_.size()) {
        D_.clear();
        for (int N = 0; N < C_left_.size(); N++) {
            std::stringstream s;
            s << "D " << N << " (SO)";
            const auto& row_sizes = C_left_[N]->rowspi();
            const auto& col_sizes = (N < C_right_.size()) ? C_right_[N]->rowspi() : row_sizes;
            D_.push_back(std::make_shared<ComplexMatrix>(s.str(), row_sizes, col_sizes));
        }
    }

    for (int N = 0; N < D_.size(); ++N) {
        D_[N]->zero();

        auto const& Cl = *C_left_[N];
        // compute() aliases C_right_ = C_left_ when only C_left is filled; allow the same here.
        auto const& Cr = (N < C_right_.size()) ? *C_right_[N] : Cl;

        if (Cl.nirrep() != Cr.nirrep()) {
            throw PSIEXCEPTION("ComplexJK::compute_D: C_left/C_right must have same number of irreps.");
        }

        for (int h = 0; h < D_[N]->nirrep(); ++h) {
            int nsol = Cl.rowdim(h);
            int nsor = Cr.rowdim(h);
            int nocc = Cl.coldim(h);
            if (nocc != Cr.coldim(h)) {
                throw PSIEXCEPTION("ComplexJK::compute_D: C_left/C_right occupied dimensions must match.");
            }
            if (!nsol || !nsor || !nocc) continue;

            auto const& Cl_h = Cl.get(h);
            auto const& Cr_h = Cr.get(h);
            auto& D_h = D_[N]->get(h);  // allocates (nsol x nsor)

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
        for (auto J : J_) J->zero();
    }
    if (do_K_) {
        for (auto K : K_) K->zero();
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

    // Figure out the symmetry and which codes will stay in C1 symmetry
    input_symmetry_cast_map_.clear();
    for (size_t i = 0; i < C_left_.size(); i++) {
        if (C_left_[i]->nirrep() != C_right_[i]->nirrep()) {
            throw PSIEXCEPTION("ComplexJK: C_left/C_right irrep mismatch!");
        }

        if ((AO2USO_->nirrep() == 1) && (C_left_[i]->nirrep() == 1)) {
            // Everything in C1, nothing to do
            input_symmetry_cast_map_.push_back(false);
        } else if (C_left_[i]->nirrep() == AO2USO_->nirrep()) {
            // We match symmetry, does this code use C1?
            if (C1()) {
                input_symmetry_cast_map_.push_back(true);
            } else {
                input_symmetry_cast_map_.push_back(false);
            }
        } else if ((C_left_[i]->nirrep() == 1) && C1()) {
            // Code uses C1, nothing to do
            input_symmetry_cast_map_.push_back(false);
        } else {
            throw PSIEXCEPTION("ComplexJK: Input orbital irrep mismatch!");
        }
    }

    // Construct the densities
    timer_on("JK: D");
    compute_D();
    timer_off("JK: D");

    bool need_cast = !input_symmetry_cast_map_.empty() &&
                     std::any_of(input_symmetry_cast_map_.begin(), input_symmetry_cast_map_.end(),
                                [](bool v) { return v; });

    if (need_cast) {
        timer_on("JK: USO2AO");
        USO2AO();
        timer_off("JK: USO2AO");
    } else {
        allocate_JK();
        // When not casting, the AO pointers alias the SO objects
        D_ao_ = D_;
        J_ao_ = J_;
        K_ao_ = K_;
    }

    timer_on("JK: JK");
    compute_JK();
    timer_off("JK: JK");

    if (need_cast) {
        timer_on("JK: AO2USO");
        AO2USO();
        timer_off("JK: AO2USO");
    }

    // Clear AO aliases so they don't hold refs
    D_ao_.clear();
    J_ao_.clear();
    K_ao_.clear();

    if (debug_ > 6) {
        for (size_t N = 0; N < C_left_.size(); N++) {
            C_left_[N]->print("outfile");
            C_right_[N]->print("outfile");
            D_[N]->print("outfile");
            if (N < J_.size()) J_[N]->print("outfile");
            if (N < K_.size()) K_[N]->print("outfile");
        }
    }

    if (lr_symmetric_) {
        C_right_.clear();
    }
}

void ComplexJK::USO2AO() {
    // Allocate D_ao_, J_ao_, K_ao_
    int nao = AO2USO_->rowspi()[0];
    int nirrep = AO2USO_->nirrep();

    D_ao_.clear();
    J_ao_.clear();
    K_ao_.clear();

    for (size_t N = 0; N < D_.size(); N++) {
        D_ao_.push_back(std::make_shared<ComplexMatrix>("D (AO)", nao, nao));
        D_ao_[N]->zero();
    }

    for (size_t N = 0; N < D_.size() && do_J_; N++) {
        J_ao_.push_back(std::make_shared<ComplexMatrix>("J (AO)", nao, nao));
    }

    for (size_t N = 0; N < D_.size() && do_K_; N++) {
        K_ao_.push_back(std::make_shared<ComplexMatrix>("K (AO)", nao, nao));
    }

    // Transform each D from SO (symmetry-blocked) to AO (C1)
    for (size_t N = 0; N < D_.size(); N++) {
        if (!input_symmetry_cast_map_[N]) {
            // Already C1 — just copy the single block
            D_ao_[N]->get(0) = D_[N]->get(0);
            continue;
        }

        // D_ao = sum_h U[h] * D_so[h] * U[h]^T
        using einsums::blas::gemm;
        for (int h = 0; h < nirrep; h++) {
            int nsol = AO2USO_->colspi()[h];
            if (nsol == 0) continue;

            if (D_[N]->rowdim(h) == 0) continue;

            const auto& D_so = D_[N]->get(h);
            const auto& U_h = AO2USO_->pointer(h);

            // temp (nsol x nao) = D_so (nsol x nsol) * U_h^T (nsol x nao)
            // Since U_h is real, we need to copy U_h^T into a complex temp
            int nsor = AO2USO_->colspi()[h];
            einsums::Tensor<std::complex<double>, 2> temp("temp", nsol, nao);
            for (int i = 0; i < nsol; i++) {
                for (int j = 0; j < nao; j++) {
                    temp(i, j) = std::complex<double>(U_h[j][i], 0.0);  // U_h is (nao x nsol), transpose
                }
            }

            // temp (nsol x nao) = D_so * temp
            einsums::Tensor<std::complex<double>, 2> temp2("temp2", nsol, nao);
            gemm('n', 'n', nsol, nao, nsor,
                 std::complex<double>{1.0}, D_so.data(), static_cast<int>(D_so.stride(0)),
                 temp.data(), static_cast<int>(temp.stride(0)),
                 std::complex<double>{0.0}, temp2.data(), static_cast<int>(temp2.stride(0)));

            // D_ao += U_h * temp2
            auto& D_ao_block = D_ao_[N]->get(0);
            einsums::Tensor<std::complex<double>, 2> U_complex("U_c", nao, nsol);
            for (int i = 0; i < nao; i++) {
                for (int j = 0; j < nsol; j++) {
                    U_complex(i, j) = std::complex<double>(U_h[i][j], 0.0);
                }
            }

            einsums::Tensor<std::complex<double>, 2> contrib("contrib", nao, nao);
            gemm('n', 'n', nao, nao, nsol,
                 std::complex<double>{1.0}, U_complex.data(), static_cast<int>(U_complex.stride(0)),
                 temp2.data(), static_cast<int>(temp2.stride(0)),
                 std::complex<double>{0.0}, contrib.data(), static_cast<int>(contrib.stride(0)));

            // Add contribution
            for (int i = 0; i < nao; i++) {
                for (int j = 0; j < nao; j++) {
                    D_ao_block(i, j) += contrib(i, j);
                }
            }
        }
    }

    // Allocate SO J/K (these will be filled by AO2USO later)
    allocate_JK();
}

void ComplexJK::AO2USO() {
    int nao = AO2USO_->rowspi()[0];
    int nirrep = AO2USO_->nirrep();

    for (size_t N = 0; N < D_.size(); N++) {
        if (!input_symmetry_cast_map_[N]) {
            // Already C1 — just copy back
            if (do_J_ && N < J_ao_.size() && N < J_.size()) {
                J_[N]->get(0) = J_ao_[N]->get(0);
            }
            if (do_K_ && N < K_ao_.size() && N < K_.size()) {
                K_[N]->get(0) = K_ao_[N]->get(0);
            }
            continue;
        }

        using einsums::blas::gemm;
        for (int h = 0; h < nirrep; h++) {
            int nsol = AO2USO_->colspi()[h];
            if (nsol == 0) continue;

            const auto& U_h = AO2USO_->pointer(h);
            int nsor = nsol;

            if (do_J_ && N < J_ao_.size() && N < J_.size()) {
                const auto& J_ao_block = J_ao_[N]->get(0);

                // temp (nao x nsor) = J_ao (nao x nao) * U_h (nao x nsor)
                einsums::Tensor<std::complex<double>, 2> U_comp("U_c", nao, nsor);
                for (int i = 0; i < nao; i++) {
                    for (int j = 0; j < nsor; j++) {
                        U_comp(i, j) = std::complex<double>(U_h[i][j], 0.0);
                    }
                }

                einsums::Tensor<std::complex<double>, 2> temp("temp", nao, nsor);
                gemm('n', 'n', nao, nsor, nao,
                     std::complex<double>{1.0}, J_ao_block.data(), static_cast<int>(J_ao_block.stride(0)),
                     U_comp.data(), static_cast<int>(U_comp.stride(0)),
                     std::complex<double>{0.0}, temp.data(), static_cast<int>(temp.stride(0)));

                // J_so[h] = U_h^T * temp
                einsums::Tensor<std::complex<double>, 2> U_comp_T("U_cT", nsol, nao);
                for (int i = 0; i < nsol; i++) {
                    for (int j = 0; j < nao; j++) {
                        U_comp_T(i, j) = std::complex<double>(U_h[j][i], 0.0);
                    }
                }

                auto& J_so_h = J_[N]->get(h);
                gemm('n', 'n', nsol, nsor, nao,
                     std::complex<double>{1.0}, U_comp_T.data(), static_cast<int>(U_comp_T.stride(0)),
                     temp.data(), static_cast<int>(temp.stride(0)),
                     std::complex<double>{0.0}, J_so_h.data(), static_cast<int>(J_so_h.stride(0)));
            }

            if (do_K_ && N < K_ao_.size() && N < K_.size()) {
                const auto& K_ao_block = K_ao_[N]->get(0);

                einsums::Tensor<std::complex<double>, 2> U_comp("U_c", nao, nsor);
                for (int i = 0; i < nao; i++) {
                    for (int j = 0; j < nsor; j++) {
                        U_comp(i, j) = std::complex<double>(U_h[i][j], 0.0);
                    }
                }

                einsums::Tensor<std::complex<double>, 2> temp("temp", nao, nsor);
                gemm('n', 'n', nao, nsor, nao,
                     std::complex<double>{1.0}, K_ao_block.data(), static_cast<int>(K_ao_block.stride(0)),
                     U_comp.data(), static_cast<int>(U_comp.stride(0)),
                     std::complex<double>{0.0}, temp.data(), static_cast<int>(temp.stride(0)));

                einsums::Tensor<std::complex<double>, 2> U_comp_T("U_cT", nsol, nao);
                for (int i = 0; i < nsol; i++) {
                    for (int j = 0; j < nao; j++) {
                        U_comp_T(i, j) = std::complex<double>(U_h[j][i], 0.0);
                    }
                }

                auto& K_so_h = K_[N]->get(h);
                gemm('n', 'n', nsol, nsor, nao,
                     std::complex<double>{1.0}, U_comp_T.data(), static_cast<int>(U_comp_T.stride(0)),
                     temp.data(), static_cast<int>(temp.stride(0)),
                     std::complex<double>{0.0}, K_so_h.data(), static_cast<int>(K_so_h.stride(0)));
            }
        }
    }
}

}  // namespace psi
