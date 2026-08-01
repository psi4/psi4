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
    if (screening_type == "DENSITY") {
        throw PSIEXCEPTION("ComplexDirectJK does not support SCREENING=DENSITY yet.");
    }
    do_csam_ = (screening_type == "CSAM");
    computed_shells_per_iter_["Quartets"] = {};
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
    const size_t nbf = static_cast<size_t>(primary_->nbf());

    for (size_t N = 0; N < D_.size(); N++) {
        if (D_[N]->grid_size(0) != 1)
            throw PSIEXCEPTION("Non-C1 symmetries not allowed with ComplexJK and SCF_TYPE DIRECT");

        if (!(do_J_ && do_K_)) {
            // TODO: figure out later
            throw PSIEXCEPTION("Both J and K must be computed with ComplexJK and SCF_TYPE DIRECT");
        }

        const auto& D_ref = D_[N]->tile(0, 0);
        const size_t dim = D_ref.dim(0);
        auto& J_out = J_[N]->tile(0, 0);
        auto& K_out = K_[N]->tile(0, 0);

        auto ints = std::shared_ptr<TwoBodyAOInt>(factory->eri());

        if (dim == nbf) {
            // Plain (non spin-blocked) complex density
            build_JK_matrices(ints, D_ref, J_out, K_out);
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
            ComplexT scratch("scratch", nbf, nbf);

            build_JK_matrices(ints, D_tot, J_tot, scratch);
            build_JK_matrices(ints, D_aa, scratch, K_aa);
            build_JK_matrices(ints, D_bb, scratch, K_bb);
            build_JK_matrices(ints, D_ab, scratch, K_ab);
            build_JK_matrices(ints, D_ba, scratch, K_ba);

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

void ComplexDirectJK::build_JK_matrices(std::shared_ptr<TwoBodyAOInt> ints, const ComplexT& D, ComplexT& J, ComplexT& K) {
    // TODO: Something like below
    // bool build_J = (!J.empty());
    // bool build_K = (!K.empty());
    timer_on("build_JK_matrices");

    J.zero();
    K.zero();

    // => Atomic Task Blocking <= //
    // One task = all shells on one atom. Shells stay in basis order (identity map).

    const int nshell = primary_->nshell();
    std::vector<int> task_shells;
    std::vector<int> task_starts;

    int atomic_ind = -1;
    for (int P = 0; P < nshell; P++) {
        if (primary_->shell(P).ncenter() > atomic_ind) {
            task_starts.push_back(P);
            atomic_ind++;
        }
        task_shells.push_back(P);
    }
    task_starts.push_back(nshell);

    const size_t ntask = task_starts.size() - 1;

    std::vector<int> task_offsets;
    task_offsets.push_back(0);
    for (int P2 = 0; P2 < nshell; P2++) {
        task_offsets.push_back(task_offsets[P2] + primary_->shell(task_shells[P2]).nfunction());
    }

    size_t max_task = 0;
    for (size_t task = 0; task < ntask; task++) {
        size_t size = 0;
        for (int P2 = task_starts[task]; P2 < task_starts[task + 1]; P2++) {
            size += primary_->shell(task_shells[P2]).nfunction();
        }
        max_task = std::max(max_task, size);
    }

    // Significant atom-task pairs (PQ|-- full, no uniqueness yet)
    std::vector<std::pair<int, int>> task_pairs;
    for (size_t Ptask = 0; Ptask < ntask; Ptask++) {
        for (size_t Qtask = 0; Qtask < ntask; Qtask++) {
            bool found = false;
            for (int P2 = task_starts[Ptask]; P2 < task_starts[Ptask + 1]; P2++) {
                for (int Q2 = task_starts[Qtask]; Q2 < task_starts[Qtask + 1]; Q2++) {
                    if (ints->shell_pair_significant(task_shells[P2], task_shells[Q2])) {
                        found = true;
                        task_pairs.emplace_back(Ptask, Qtask);
                        break;
                    }
                }
                if (found) break;
            }
        }
    }
    const size_t ntask_pair = task_pairs.size();

    // Local AO blocks for the current atom quartet: J_PQ and K_PR
    ComplexT JT("JT", max_task, max_task);
    ComplexT KT("KT", max_task, max_task);

    size_t computed_shells = 0L;

    for (size_t task_pq = 0; task_pq < ntask_pair; task_pq++) {
        for (size_t task_rs = 0; task_rs < ntask_pair; task_rs++) {
            const int Ptask = task_pairs[task_pq].first;
            const int Qtask = task_pairs[task_pq].second;
            const int Rtask = task_pairs[task_rs].first;
            const int Stask = task_pairs[task_rs].second;

            const int P2start = task_starts[Ptask];
            const int Q2start = task_starts[Qtask];
            const int R2start = task_starts[Rtask];
            const int S2start = task_starts[Stask];

            const int nPtask = task_starts[Ptask + 1] - P2start;
            const int nQtask = task_starts[Qtask + 1] - Q2start;
            const int nRtask = task_starts[Rtask + 1] - R2start;
            const int nStask = task_starts[Stask + 1] - S2start;

            JT.zero();
            KT.zero();
            bool touched = false;

            for (int P2 = P2start; P2 < P2start + nPtask; P2++) {
                for (int Q2 = Q2start; Q2 < Q2start + nQtask; Q2++) {
                    const int P = task_shells[P2];
                    const int Q = task_shells[Q2];
                    if (!ints->shell_pair_significant(P, Q)) continue;

                    for (int R2 = R2start; R2 < R2start + nRtask; R2++) {
                        for (int S2 = S2start; S2 < S2start + nStask; S2++) {
                            const int R = task_shells[R2];
                            const int S = task_shells[S2];
                            if (!ints->shell_pair_significant(R, S)) continue;
                            if (!ints->shell_significant(P, Q, R, S)) continue;

                            if (ints->compute_shell(P, Q, R, S) == 0) continue;
                            const double* buf = ints->buffer();
                            computed_shells++;

                            const int Psize = primary_->shell(P).nfunction();
                            const int Qsize = primary_->shell(Q).nfunction();
                            const int Rsize = primary_->shell(R).nfunction();
                            const int Ssize = primary_->shell(S).nfunction();
                            const int Poff = primary_->shell(P).function_index();
                            const int Qoff = primary_->shell(Q).function_index();
                            const int Roff = primary_->shell(R).function_index();
                            const int Soff = primary_->shell(S).function_index();
                            const int Poff2 = task_offsets[P2] - task_offsets[P2start];
                            const int Qoff2 = task_offsets[Q2] - task_offsets[Q2start];
                            const int Roff2 = task_offsets[R2] - task_offsets[R2start];

                            touched = true;
                            for (int p = 0; p < Psize; p++) {
                                for (int q = 0; q < Qsize; q++) {
                                    for (int r = 0; r < Rsize; r++) {
                                        for (int s = 0; s < Ssize; s++) {
                                            const double I = *buf++;
                                            JT(p + Poff2, q + Qoff2) += D(r + Roff, s + Soff) * I;
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
            for (int P2 = 0; P2 < nPtask; P2++) {
                for (int Q2 = 0; Q2 < nQtask; Q2++) {
                    const int P = task_shells[P2start + P2];
                    const int Q = task_shells[Q2start + Q2];
                    const int Psize = primary_->shell(P).nfunction();
                    const int Qsize = primary_->shell(Q).nfunction();
                    const int Poff = primary_->shell(P).function_index();
                    const int Qoff = primary_->shell(Q).function_index();
                    const int Poff2 = task_offsets[P2start + P2] - task_offsets[P2start];
                    const int Qoff2 = task_offsets[Q2start + Q2] - task_offsets[Q2start];
                    for (int p = 0; p < Psize; p++) {
                        for (int q = 0; q < Qsize; q++) {
                            J(p + Poff, q + Qoff) += JT(p + Poff2, q + Qoff2);
                        }
                    }
                }
            }

            for (int P2 = 0; P2 < nPtask; P2++) {
                for (int R2 = 0; R2 < nRtask; R2++) {
                    const int P = task_shells[P2start + P2];
                    const int R = task_shells[R2start + R2];
                    const int Psize = primary_->shell(P).nfunction();
                    const int Rsize = primary_->shell(R).nfunction();
                    const int Poff = primary_->shell(P).function_index();
                    const int Roff = primary_->shell(R).function_index();
                    const int Poff2 = task_offsets[P2start + P2] - task_offsets[P2start];
                    const int Roff2 = task_offsets[R2start + R2] - task_offsets[R2start];
                    for (int p = 0; p < Psize; p++) {
                        for (int r = 0; r < Rsize; r++) {
                            K(p + Poff, r + Roff) += KT(p + Poff2, r + Roff2);
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
