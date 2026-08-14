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
#include "psi4/liboptions/liboptions.h"

#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libqt/qt.h"

#include <complex>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <ranges>

namespace psi {

ComplexDirectJK::ComplexDirectJK(std::shared_ptr<BasisSet> primary, Options& options)
    : ComplexJK(primary), options_(options) {
    common_init();
}

ComplexDirectJK::~ComplexDirectJK() {}

void ComplexDirectJK::common_init() {
    if (options_.get_int("INCFOCK_FULL_FOCK_EVERY") <= 0) {
        throw PSIEXCEPTION("Invalid input for option INCFOCK_FULL_FOCK_EVERY (<= 0)");
    }

    auto screening_type = options_.get_str("SCREENING");
    // if (screening_type == "DENSITY") {
    //     throw PSIEXCEPTION("ComplexDirectJK does not support SCREENING=DENSITY yet.");
    // }
    do_csam_ = (screening_type == "CSAM");
    computed_shells_per_iter_["Quartets"] = {};

    // Set up AO2USO transform 
    ComplexJK::common_init();
}

size_t ComplexDirectJK::num_computed_shells() { return num_computed_shells_; }

size_t ComplexDirectJK::memory_estimate() {
    return 0; // Effectively?
}

void ComplexDirectJK::preiterations() { /*no-op*/ }

void ComplexDirectJK::compute_JK() {
    if (!do_J_ && !do_K_) {
        outfile->Printf("\n  WARNING: ComplexDirectJK tried to compute JK, but "
                        "do_J_ and do_K_ were both set to false.\n");
        return;
    }

    num_computed_shells_ = 0L;

    auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    const int nbf = primary_->nbf();

    auto ints = std::shared_ptr<TwoBodyAOInt>(factory->eri());
    if (options_.get_str("SCREENING") == "DENSITY") ints->update_density_complex(D_ao_);

    for (size_t N = 0; N < D_.size(); N++) {
        if (!(do_J_ && do_K_)) {
            // TODO: figure out later
            throw PSIEXCEPTION("Both J and K must be computed with ComplexJK and SCF_TYPE DIRECT");
        }

        const auto& D_ref = D_ao_[N]->get(0);
        const int dim = D_ref.dim(0);
        auto& J_out = J_ao_[N]->get(0);
        auto& K_out = K_ao_[N]->get(0);

        if (dim == nbf) {
            // Plain (non spin-blocked) complex density
            build_JK_matrices<true, true>(ints, D_ref, &J_out, &K_out);
        } else if (dim == 2 * nbf) {
            // Generalized (CGHF) spin-blocked density. D and J/K are 2x2 block
            // matrices of nbf x nbf blocks: [[D_aa, D_ab], [D_ba, D_bb]].
            // A single ERI pass contracts all four spin blocks (and the summed
            // Coulomb density) so the integrals are only computed once.
            build_JK_matrices<true, true, true>(ints, D_ref, &J_out, &K_out);
        } else {
            throw PSIEXCEPTION("ComplexDirectJK: density block dimension (" + std::to_string(dim) +
                               ") must equal the number of AO basis functions (" + std::to_string(nbf) +
                               ") or twice that (GHF spin-blocked density).");
        }
    }

    if (get_bench()) {
        computed_shells_per_iter_["Quartets"].push_back(num_computed_shells());
    }
}

namespace {

// Helper method to remove the indexing madness. Produces iterable list of iterators.
// Converts Something like [0,1,4,5] into [[0, 1), [1, 4), [4, 5)].
auto partition(std::span<const std::size_t> atom_to_shell) {
return std::views::iota(std::size_t{0}, atom_to_shell.size() - 1)
     | std::views::transform([atom_to_shell](std::size_t atom) {
           return std::views::iota(atom_to_shell[atom], atom_to_shell[atom + 1]);
       });
}

// Like above but enumerated. Converts Something like [0,1,4,5] into
//   [{0, [0, 1)}, {1, [1, 4)}, {2, [4, 5)}].
auto partition_with_idx(std::span<const std::size_t> atom_to_shell) {
return std::views::iota(std::size_t{0}, atom_to_shell.size() - 1)
     | std::views::transform([atom_to_shell](std::size_t atom) {
           return std::pair{atom, std::views::iota(atom_to_shell[atom], atom_to_shell[atom + 1])};
       });
}

}

// The template allows us to provide the most optimized method at compile time
// using constexpr. No longer using scratch matrices. Benchmarking shows this
// being only marginally faster (<5%), but the code smells less.
//
// SpinBlocked selects between a plain complex density (nbf x nbf) and a
// generalized (CGHF) spin-blocked density (2 nbf x 2 nbf). It is a compile-time
// switch: both paths share the atom/shell-blocked skeleton and differ only in
// how the density is decomposed and how J/K are accumulated and striped back,
// with no runtime branching inside the hot loops.
template<bool FillJ, bool FillK, bool SpinBlocked>
void ComplexDirectJK::build_JK_matrices(std::shared_ptr<TwoBodyAOInt> ints, const ComplexT& D, ComplexT* J, ComplexT* K) {
    timer_on("build_JK_matrices");

    [[maybe_unused]] const int nbf = primary_->nbf();

    // Spin-blocked (CGHF) density decomposition. A generalized density is
    // D = [[D_aa, D_ab], [D_ba, D_bb]]; the Coulomb term only needs the summed
    // charge density D_aa + D_bb, while exchange needs each block individually.
    // The blocks are materialized as contiguous nbf x nbf copies because the
    // corresponding blocks of the 2 nbf x 2 nbf density are strided.
    [[maybe_unused]] std::optional<ComplexT> D_aa, D_bb, D_ab, D_ba, D_tot;
    if constexpr (SpinBlocked) {
        D_aa.emplace("D_aa", nbf, nbf);
        D_bb.emplace("D_bb", nbf, nbf);
        D_ab.emplace("D_ab", nbf, nbf);
        D_ba.emplace("D_ba", nbf, nbf);
        D_tot.emplace("D_tot", nbf, nbf);
        for (size_t p = 0; p < nbf; p++) {
            for (size_t q = 0; q < nbf; q++) {
                (*D_aa)(p, q) = D(p, q);
                (*D_bb)(p, q) = D(p + nbf, q + nbf);
                (*D_ab)(p, q) = D(p, q + nbf);
                (*D_ba)(p, q) = D(p + nbf, q);
                (*D_tot)(p, q) = (*D_aa)(p, q) + (*D_bb)(p, q);
            }
        }
    }

    if constexpr(FillJ) J->zero();
    if constexpr(FillK) K->zero();

    // => Atomic Task Blocking <= //
    // One task = all shells on one center. Shells stay in basis order (identity map).

    // shell index that starts each center
    std::vector<size_t> atom_to_shell;
    // boundaries of each partition
    std::vector<size_t> shell_to_basis;

    /*******************************************************
     * Consider H2O with STO-3G...                         *
     *                                                     *
     *   3 atoms:       H  ,         O           ,   H     *
     *   5 shells:    [ S ], [ S , S ,    P     ], [ S ]   *
     *   7 functions: [[s]], [[s],[s],[px,py,pz]], [[s]]   *
     *                                                     *
     * atom_to_shell  [ 0,           1,              4, 5] *
     * shell_to_basis [ 0,     1,  2,     3,         6, 7] *
     *                                                     *
     *******************************************************/

    // => Helper functions that need internal variables <= //

    auto nfunctions_in_shell = [this](const size_t& shell_idx) {
        return this->primary_->shell(shell_idx).nfunction();
    };

    auto function_index_of_shell = [this](const size_t& shell_idx) {
        return this->primary_->shell(shell_idx).function_index();
    };

    // Returns the basis function index that begins shell ``shell_idx`` where 0
    // indicates the first basis function of the first shell of atom ``task``.
    auto basis_index_of_shell_from_atom = [&atom_to_shell, &shell_to_basis](const size_t& shell_idx, const size_t& task) {
        return shell_to_basis[shell_idx] - shell_to_basis[atom_to_shell[task]];
    };

    // Returns a view of shell indices for given center idx.
    auto shells_on_center = [&atom_to_shell](size_t task) {
        return std::views::iota(atom_to_shell[task], atom_to_shell[task+1]);
    };

    // => Welcome to the jungle <=

    {
        // Total number of shells for all centers
        const int nshell = primary_->nshell();

        int atomic_ind = -1;

        size_t total_nfuncs = 0;
        shell_to_basis.push_back(0);
        for (int P = 0; P < nshell; P++) {
            const auto& shell = primary_->shell(P);

            total_nfuncs += shell.nfunction();
            shell_to_basis.push_back(total_nfuncs);

            if (shell.ncenter() > atomic_ind) {
                atom_to_shell.push_back(P);
                atomic_ind++;
            }
        }
        atom_to_shell.push_back(nshell);
    }

    // Largest number of functions for any center. Used to allocate temporary
    // matrix for J and K.
    size_t max_nfuncs_per_center = 0;
    for (auto centers : partition(atom_to_shell)) {
        size_t size = 0;
        for (size_t shell : centers) {
            size += primary_->shell(shell).nfunction();
        }
        max_nfuncs_per_center = std::max(max_nfuncs_per_center, size);
    }

    // Local AO blocks for the current atom quartet: J_PQ and K_PR
    ComplexT JT("JT", max_nfuncs_per_center, max_nfuncs_per_center);
    // Exchange accumulators: a plain density contracts into a single block, a
    // spin-blocked density into one block per spin component (aa/ab/ba/bb).
    [[maybe_unused]] ComplexT KT("KT", max_nfuncs_per_center, max_nfuncs_per_center);
    [[maybe_unused]] ComplexT KaaT("KaaT", max_nfuncs_per_center, max_nfuncs_per_center);
    [[maybe_unused]] ComplexT KabT("KabT", max_nfuncs_per_center, max_nfuncs_per_center);
    [[maybe_unused]] ComplexT KbaT("KbaT", max_nfuncs_per_center, max_nfuncs_per_center);
    [[maybe_unused]] ComplexT KbbT("KbbT", max_nfuncs_per_center, max_nfuncs_per_center);

    // Significant atom/task pairs (PQ|-- full, no uniqueness yet)
    std::vector<std::pair<size_t, size_t>> significant_pairs;
    for (auto [Patom, Pshells] : partition_with_idx(atom_to_shell)) {
        for (auto [Qatom, Qshells] : partition_with_idx(atom_to_shell)) {
            bool found_significant_pair = false;
            for (size_t Ps : Pshells) {
                for (size_t Qs : Qshells) {
                    if (ints->shell_pair_significant(Ps, Qs)) {
                        found_significant_pair = true;
                        significant_pairs.emplace_back(Patom, Qatom);
                        break;
                    }
                }
                if (found_significant_pair) break;
            }
        }
    }

    size_t computed_shells = 0L;
    for (const auto [Patom, Qatom] : significant_pairs) {
        for (const auto [Ratom, Satom] : significant_pairs) {
            if constexpr(FillJ) JT.zero();
            if constexpr(FillK) {
                if constexpr (SpinBlocked) {
                    KaaT.zero();
                    KabT.zero();
                    KbaT.zero();
                    KbbT.zero();
                } else {
                    KT.zero();
                }
            }

            bool touched = false;

            // Process the shells
            for (auto Ps : shells_on_center(Patom) ) {
                for (auto Qs : shells_on_center(Qatom) ) {
                    if (!ints->shell_pair_significant(Ps, Qs)) continue;
                    for (auto Rs : shells_on_center(Ratom)) {
                        for (auto Ss : shells_on_center(Satom)) {
                            if (!ints->shell_pair_significant(Rs, Ss)) continue;
                            if (!ints->shell_significant(Ps, Qs, Rs, Ss)) continue;

                            // Compute ERI
                            if (ints->compute_shell(Ps, Qs, Rs, Ss) == 0) continue;
                            const double* buf = ints->buffer();

                            // Don't forget to like, share, subscribe and compute that shell!
                            computed_shells++;

                            const size_t Psize = nfunctions_in_shell(Ps);
                            const size_t Qsize = nfunctions_in_shell(Qs);
                            const size_t Rsize = nfunctions_in_shell(Rs);
                            const size_t Ssize = nfunctions_in_shell(Ss);

                            // Global basis function index. 0 means the first function
                            // of the first shell of the first atom. aka AO index.
                            const size_t Qao = function_index_of_shell(Qs);
                            const size_t Rao = function_index_of_shell(Rs);
                            const size_t Sao = function_index_of_shell(Ss);

                            // basis index starting shell Ps from atom Patom.
                            // Let's call these "lo" for local orbital. Like
                            // AO but where 0 means the first AO of the given atom.
                            const size_t Plo = basis_index_of_shell_from_atom(Ps, Patom);
                            const size_t Qlo = basis_index_of_shell_from_atom(Qs, Qatom);
                            const size_t Rlo = basis_index_of_shell_from_atom(Rs, Ratom);

                            touched = true;
                            for (size_t p = 0; p < Psize; p++) {
                                for (size_t q = 0; q < Qsize; q++) {
                                    for (size_t r = 0; r < Rsize; r++) {
                                        for (size_t s = 0; s < Ssize; s++) {
                                            const double I = *buf++;
                                            if constexpr(FillJ) {
                                                if constexpr (SpinBlocked)
                                                    JT(p + Plo, q + Qlo) += (*D_tot)(r + Rao, s + Sao) * I;
                                                else
                                                    JT(p + Plo, q + Qlo) += D(r + Rao, s + Sao) * I;
                                            }
                                            if constexpr(FillK) {
                                                if constexpr (SpinBlocked) {
                                                    KaaT(p + Plo, r + Rlo) += (*D_aa)(q + Qao, s + Sao) * I;
                                                    KabT(p + Plo, r + Rlo) += (*D_ab)(q + Qao, s + Sao) * I;
                                                    KbaT(p + Plo, r + Rlo) += (*D_ba)(q + Qao, s + Sao) * I;
                                                    KbbT(p + Plo, r + Rlo) += (*D_bb)(q + Qao, s + Sao) * I;
                                                } else {
                                                    KT(p + Plo, r + Rlo) += D(q + Qao, s + Sao) * I;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (!touched) continue;

            // Stripe task-local J_PQ / K_PR into global matrices
            if constexpr(FillJ) {
                for (auto Ps : shells_on_center(Patom) ) {
                    for (auto Qs : shells_on_center(Qatom) ) {
                        const size_t Psize = nfunctions_in_shell(Ps);
                        const size_t Qsize = nfunctions_in_shell(Qs);

                        const size_t Pao = function_index_of_shell(Ps);
                        const size_t Qao = function_index_of_shell(Qs);

                        const size_t Plo = basis_index_of_shell_from_atom(Ps, Patom);
                        const size_t Qlo = basis_index_of_shell_from_atom(Qs, Qatom);
                        for (int p = 0; p < Psize; p++) {
                            for (int q = 0; q < Qsize; q++) {
                                if constexpr (SpinBlocked) {
                                    // J_aa and J_bb are identical: J only depends
                                    // on the summed charge density D_aa + D_bb.
                                    (*J)(p + Pao, q + Qao) += JT(p + Plo, q + Qlo);
                                    (*J)(p + Pao + nbf, q + Qao + nbf) += JT(p + Plo, q + Qlo);
                                } else {
                                    (*J)(p + Pao, q + Qao) += JT(p + Plo, q + Qlo);
                                }
                            }
                        }
                    }
                }
            }

            if constexpr(FillK) {
                for (auto Ps : shells_on_center(Patom) ) {
                    for (auto Rs : shells_on_center(Ratom) ) {
                        const size_t Psize = nfunctions_in_shell(Ps);
                        const size_t Rsize = nfunctions_in_shell(Rs);

                        const size_t Pao = function_index_of_shell(Ps);
                        const size_t Rao = function_index_of_shell(Rs);

                        const size_t Plo = basis_index_of_shell_from_atom(Ps, Patom);
                        const size_t Rlo = basis_index_of_shell_from_atom(Rs, Ratom);
                        for (int p = 0; p < Psize; p++) {
                            for (int r = 0; r < Rsize; r++) {
                                if constexpr (SpinBlocked) {
                                    (*K)(p + Pao, r + Rao) += KaaT(p + Plo, r + Rlo);
                                    (*K)(p + Pao, r + Rao + nbf) += KabT(p + Plo, r + Rlo);
                                    (*K)(p + Pao + nbf, r + Rao) += KbaT(p + Plo, r + Rlo);
                                    (*K)(p + Pao + nbf, r + Rao + nbf) += KbbT(p + Plo, r + Rlo);
                                } else {
                                    (*K)(p + Pao, r + Rao) += KT(p + Plo, r + Rlo);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    num_computed_shells_ += computed_shells;
    timer_off("build_JK_matrices");
}

void ComplexDirectJK::postiterations() { /*no-op*/ }

void ComplexDirectJK::print_header() const {
    std::string screen_type = options_.get_str("SCREENING");
    if (print_) {
        outfile->Printf("  ==> ComplexDirectJK: Integral-Direct J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    Memory [MiB]:      %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf("    Screening Type:    %11s\n", screen_type.c_str());
        outfile->Printf("    Screening Cutoff:  %11.0E\n", cutoff_);
        outfile->Printf("\n");
    }
}

} // namespace psi
