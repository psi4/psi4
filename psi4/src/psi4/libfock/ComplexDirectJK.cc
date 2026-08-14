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
            // Generalized (CGHF) spin-blocked density. (D and J/K) are 2x2 block
            // matrices of nbf x nbf blocks: [[D_aa, D_ab], [D_ba, D_bb]]

            // This still leaves room for improvement in terms of uniqueness.
            // D_ab is often adjoint of D_ba, for example. D_aa and D_bb are
            // often Hermitian themselves.

            // Splitting up into separate tensors is preferred. These are contiguous
            // where Views would not be in general.
            ComplexT D_aa("D_aa", nbf, nbf), D_bb("D_bb", nbf, nbf);
            ComplexT D_ab("D_ab", nbf, nbf), D_ba("D_ba", nbf, nbf);
            for (size_t p = 0; p < nbf; p++) {
                for (size_t q = 0; q < nbf; q++) {
                    D_aa(p, q) = D_ref(p, q);
                    D_bb(p, q) = D_ref(p + nbf, q + nbf);
                    D_ab(p, q) = D_ref(p, q + nbf);
                    D_ba(p, q) = D_ref(p + nbf, q);
                }
            }

            ComplexT D_tot = D_aa;
            for (size_t p = 0; p < nbf; p++)
                for (size_t q = 0; q < nbf; q++) D_tot(p, q) += D_bb(p, q);

            ComplexT J_tot("J_tot", nbf, nbf);
            ComplexT K_aa("K_aa", nbf, nbf), K_bb("K_bb", nbf, nbf);
            ComplexT K_ab("K_ab", nbf, nbf), K_ba("K_ba", nbf, nbf);

            build_JK_matrices<true, false>(ints, D_tot, &J_tot, nullptr);
            build_JK_matrices<false, true>(ints, D_aa, nullptr, &K_aa);
            build_JK_matrices<false, true>(ints, D_bb, nullptr, &K_bb);
            build_JK_matrices<false, true>(ints, D_ab, nullptr, &K_ab);
            build_JK_matrices<false, true>(ints, D_ba, nullptr, &K_ba);

            J_out.zero();
            K_out.zero();
            for (size_t p = 0; p < nbf; p++) {
                for (size_t q = 0; q < nbf; q++) {
                    J_out(p, q) = J_tot(p, q);
                    J_out(p + nbf, q + nbf) = J_tot(p, q);
                    K_out(p, q) = K_aa(p, q);
                    K_out(p + nbf, q + nbf) = K_bb(p, q);
                    K_out(p, q + nbf) = K_ab(p, q);
                    K_out(p + nbf, q) = K_ba(p, q);
                }
            }
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

// Helper method to remove the indexing maddness. Produces iterable list of iterators.
// Converts Something like starts = [0,1,4,5] into [[0, 1), [1, 4), [4, 5)].
auto partition(std::span<const std::size_t> starts) {
return std::views::iota(std::size_t{0}, starts.size() - 1)
     | std::views::transform([starts](std::size_t i) {
           return std::views::iota(starts[i], starts[i + 1]);
       });
}

auto partition_with_idx(std::span<const std::size_t> starts) {
return std::views::iota(std::size_t{0}, starts.size() - 1)
     | std::views::transform([starts](std::size_t i) {
           return std::pair{i, std::views::iota(starts[i], starts[i + 1])};
       });
}

}

// The template allows for us to provide the most optimized method at compile time
// using constexpr. No longer using scratch matrices. Benchmarking shows this
// being only marginally faster (<5%), but the code smells less.
template<bool FillJ, bool FillK>
void ComplexDirectJK::build_JK_matrices(std::shared_ptr<TwoBodyAOInt> ints, const ComplexT& D, ComplexT* J, ComplexT* K) {
    timer_on("build_JK_matrices");

    if constexpr(FillJ) J->zero();
    if constexpr(FillK) K->zero();

    // => Helpers <= //

    auto nfunctions_in_shell = [this](const size_t& shell_idx) {
        return this->primary_->shell(shell_idx).nfunction();
    };

    auto function_index_of_shell = [this](const size_t& shell_idx) {
        return this->primary_->shell(shell_idx).function_index();
    };


    // => Atomic Task Blocking <= //
    // One task = all shells on one center. Shells stay in basis order (identity map).

    // shell index that starts each center
    std::vector<size_t> center_starts;
    // boundaries of each partition
    std::vector<size_t> center_offsets;

    {
        // Total number of shells for all centers
        const int nshell = primary_->nshell();

        int atomic_ind = -1;

        size_t total_nfuncs = 0;
        center_offsets.push_back(0);
        for (int P = 0; P < nshell; P++) {
            const auto& shell = primary_->shell(P);

            total_nfuncs += shell.nfunction();
            center_offsets.push_back(total_nfuncs);

            if (shell.ncenter() > atomic_ind) {
                center_starts.push_back(P);
                atomic_ind++;
            }
        }
        center_starts.push_back(nshell);
    }

    // Number of centers (atoms and ghosts)
    const size_t ncenter = center_starts.size() - 1;

    // largest number of functions for any center.
    size_t max_nfuncs_per_center = 0;
    for (auto centers : partition(center_starts)) {
        size_t size = 0;
        for (size_t shell : centers) {
            size += primary_->shell(shell).nfunction();
        }
        max_nfuncs_per_center = std::max(max_nfuncs_per_center, size);
    }

    // Significant atom-task pairs (PQ|-- full, no uniqueness yet)
    std::vector<std::pair<size_t, size_t>> significant_pairs;
    for (auto [P_center_idx, P_shells] : partition_with_idx(center_starts)) {
        for (auto [Q_center_idx, Q_shells] : partition_with_idx(center_starts)) {
            bool found_significant_pair = false;
            for (size_t P_shell : P_shells) {
                for (size_t Q_shell : Q_shells) {
                    if (ints->shell_pair_significant(P_shell, Q_shell)) {
                        found_significant_pair = true;
                        significant_pairs.emplace_back(P_center_idx, Q_center_idx);
                        break;
                    }
                }
                if (found_significant_pair) break;
            }
        }
    }

    // Local AO blocks for the current atom quartet: J_PQ and K_PR
    ComplexT JT("JT", max_nfuncs_per_center, max_nfuncs_per_center);
    ComplexT KT("KT", max_nfuncs_per_center, max_nfuncs_per_center);

    size_t computed_shells = 0L;

    for (const auto [Ptask, Qtask] : significant_pairs) {
        for (const auto [Rtask, Stask] : significant_pairs) {
            // auto shells_on_center = [&](size_t task) {
            //     return std::views::iota(center_starts[task], center_starts[task+1])
            //         | std::views::transform([&](size_t shell_idx){
            //             return std::pair{shell_idx, center_offsets[shell_idx] - center_offsets[center_starts[task]]};
            //     });
            // };
            // Returns a view of shell indices for given center idx.
            auto shells_on_center = [&](size_t task) {
                return std::views::iota(center_starts[task], center_starts[task+1]);
            };

            // Global shell index for Ptask. 0 means the first shell of the first atom.
            const int P2start = center_starts[Ptask];
            const int Q2start = center_starts[Qtask];
            const int R2start = center_starts[Rtask];
            const int S2start = center_starts[Stask];

            // Number of shells that will be processed for this atom (Ptask)
            const int nPtask = center_starts[Ptask + 1] - P2start;
            const int nQtask = center_starts[Qtask + 1] - Q2start;
            const int nRtask = center_starts[Rtask + 1] - R2start;
            const int nStask = center_starts[Stask + 1] - S2start;

            if constexpr(FillJ) JT.zero();
            if constexpr(FillK) KT.zero();

            bool touched = false;

            // Process the shells
            // for (auto [P2, Poff2] : shells_on_center(Ptask) ) {
            //     for (auto [Q2, Qoff2] : shells_on_center(Qtask) ) {
            //         if (!ints->shell_pair_significant(P2, Q2)) continue;
            //
            //         for (auto [R2, Roff2] : shells_on_center(Rtask)) {
            //             for (auto [S2, Soff2] : shells_on_center(Stask)) {

            for (auto P2 : shells_on_center(Ptask) ) {
                for (auto Q2 : shells_on_center(Qtask) ) {
                    if (!ints->shell_pair_significant(P2, Q2)) continue;

                    for (auto R2 : shells_on_center(Rtask)) {
                        for (auto S2 : shells_on_center(Stask)) {

                            if (!ints->shell_pair_significant(R2, S2)) continue;
                            if (!ints->shell_significant(P2, Q2, R2, S2)) continue;

                            // Compute ERI
                            if (ints->compute_shell(P2, Q2, R2, S2) == 0) continue;

                            const double* buf = ints->buffer();

                            // Don't forget to like, share, subscribe and increment that variable
                            computed_shells++;

                            const size_t Psize = nfunctions_in_shell(P2);
                            const size_t Qsize = nfunctions_in_shell(Q2);
                            const size_t Rsize = nfunctions_in_shell(R2);
                            const size_t Ssize = nfunctions_in_shell(S2);

                            // Global basis function index. 0 means the first function
                            // of the first shell of the first atom. aka AO index.
                            const size_t Qoff = function_index_of_shell(Q2);
                            const size_t Roff = function_index_of_shell(R2);
                            const size_t Soff = function_index_of_shell(S2);

                            // Number of basis functions since P2start
                            const int Poff2 = center_offsets[P2] - center_offsets[P2start];
                            const int Qoff2 = center_offsets[Q2] - center_offsets[Q2start];
                            const int Roff2 = center_offsets[R2] - center_offsets[R2start];

                            touched = true;
                            for (size_t p = 0; p < Psize; p++) {
                                for (size_t q = 0; q < Qsize; q++) {
                                    for (size_t r = 0; r < Rsize; r++) {
                                        for (size_t s = 0; s < Ssize; s++) {
                                            const double I = *buf++;
                                            if constexpr(FillJ)
                                                JT(p + Poff2, q + Qoff2) += D(r + Roff, s + Soff) * I;
                                            if constexpr(FillK)
                                                KT(p + Poff2, r + Roff2) += D(q + Qoff, s + Soff) * I;
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
                for (int P2 = 0; P2 < nPtask; P2++) {
                    for (int Q2 = 0; Q2 < nQtask; Q2++) {
                        const int P = P2start + P2;
                        const int Q = Q2start + Q2;
                        const int Psize = primary_->shell(P).nfunction();
                        const int Qsize = primary_->shell(Q).nfunction();
                        const int Poff = primary_->shell(P).function_index();
                        const int Qoff = primary_->shell(Q).function_index();
                        const int Poff2 = center_offsets[P2start + P2] - center_offsets[P2start];
                        const int Qoff2 = center_offsets[Q2start + Q2] - center_offsets[Q2start];
                        for (int p = 0; p < Psize; p++) {
                            for (int q = 0; q < Qsize; q++) {
                                (*J)(p + Poff, q + Qoff) += JT(p + Poff2, q + Qoff2);
                            }
                        }
                    }
                }
            }

            if constexpr(FillK) {
                for (int P2 = 0; P2 < nPtask; P2++) {
                    for (int R2 = 0; R2 < nRtask; R2++) {
                        const int P = P2start + P2;
                        const int R = R2start + R2;
                        const int Psize = primary_->shell(P).nfunction();
                        const int Rsize = primary_->shell(R).nfunction();
                        const int Poff = primary_->shell(P).function_index();
                        const int Roff = primary_->shell(R).function_index();
                        const int Poff2 = center_offsets[P2start + P2] - center_offsets[P2start];
                        const int Roff2 = center_offsets[R2start + R2] - center_offsets[R2start];
                        for (int p = 0; p < Psize; p++) {
                            for (int r = 0; r < Rsize; r++) {
                                (*K)(p + Poff, r + Roff) += KT(p + Poff2, r + Roff2);
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
